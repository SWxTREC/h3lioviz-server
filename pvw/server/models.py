import json
import math
import os
import tempfile
from pathlib import Path

import numpy as np
import paraview.simple as pvs
from vtkmodules.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline
from vtkmodules.vtkIONetCDF import vtkNetCDFCFReader

# The dimensions of the Enlil variables, in the order the reader needs them
ENLIL_DIMENSIONS = "(longitude, latitude, radius)"

# ParaView .series file listing a run's time files and the time of each:
#   {"file-series-version": "1.0",
#    "files": [{"name": "pv-tim.0000.nc", "time": 1759644046}, ...]}
# Names are relative to the .series file and in name order, times in seconds
# since 1970 and strictly increasing. Given this file, the reader takes the
# times from it instead of opening every file to find its time (under
# ParaView 6.2, about 8 MB read from each file, twice per load). Written by
# the server on the first load of a run, or by the ingest step.
SERIES_NAME = "pv-tim.series"


def read_series(path, names_only=False):
    """
    Read a .series file.

    Returns
    -------
    list(str) or None
        The file names it lists, or None if it is missing, malformed, or
        its times are not numbers that strictly increase. With names_only,
        a file of names without times (what write_series writes when it is
        given no times) is valid too: the reader then asks each file for
        its times. Only the server's private series are written that way.
    """
    try:
        series = json.loads(Path(path).read_text())
        if series["file-series-version"] != "1.0":
            return None
        files = series["files"]
        if (
            names_only
            and isinstance(files, list)
            and files
            and all(isinstance(f, str) for f in files)
        ):
            return files
        names = [entry["name"] for entry in files]
        times = [entry["time"] for entry in series["files"]]
    # RecursionError: JSON nested too deeply for the parser
    except (OSError, ValueError, KeyError, TypeError, RecursionError):
        return None
    if not all(isinstance(name, str) for name in names):
        return None
    # An int too large for a float would make math.isfinite raise, so ints
    # are range-checked instead (a NetCDF time fits in a signed 64-bit int)
    if not all(
        not isinstance(t, bool)
        and (
            -(2**63) <= t < 2**63
            if isinstance(t, int)
            else isinstance(t, float) and math.isfinite(t)
        )
        for t in times
    ):
        return None
    if any(later <= earlier for earlier, later in zip(times, times[1:])):
        return None
    return names


def write_series(path, names, times=None):
    """
    Write a .series file, atomically: a temp file in the same directory is
    renamed over it, so a reader never sees part of one, and two sessions
    writing the same run at once each leave a whole file.

    Parameters
    ----------
    path : Path
        The .series file to write
    names : list(str)
        The file names, relative to the .series file or absolute
    times : list(float), optional
        The time of each file. Without them the reader asks each file
        for its times itself.
    """
    if times is None:
        files = list(names)
    else:
        files = [{"name": name, "time": t} for name, t in zip(names, times)]
    text = json.dumps({"file-series-version": "1.0", "files": files})
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(fd, "w") as f:
            f.write(text)
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def file_times(files):
    """
    The time of each file, as ParaView's NetCDF reader reports it.

    Returns
    -------
    list(float) or None
        One time per file, or None if a file does not hold exactly one time
        or the times do not increase in file order
    """
    reader = vtkNetCDFCFReader()
    reader.SetDimensions(ENLIL_DIMENSIONS)
    times = []
    for f in files:
        reader.SetFileName(str(f))
        reader.UpdateInformation()
        info = reader.GetOutputInformation(0)
        if not info.Has(vtkStreamingDemandDrivenPipeline.TIME_STEPS()):
            return None
        steps = info.Get(vtkStreamingDemandDrivenPipeline.TIME_STEPS())
        if steps is None or len(steps) != 1:
            return None
        times.append(steps[0])
    if any(later <= earlier for earlier, later in zip(times, times[1:])):
        return None
    return times


def private_series_path(run_dir, create=True):
    """
    Where to write a run's .series file when it cannot go in the run
    directory (not writable, or not a run of pv-tim files with one time
    each): a directory of this user's under the system temp directory.
    With create False, None if that directory does not exist yet, so only
    looking for a saved series never needs a writable temp directory.
    """
    private_dir = Path(tempfile.gettempdir()) / f"h3lioviz-series-{os.getuid()}"
    if create:
        private_dir.mkdir(mode=0o700, exist_ok=True)
    elif not private_dir.exists():
        return None
    # The name is predictable in a shared directory: refuse one that another
    # user made, or a symlink, rather than read or write series there
    if private_dir.is_symlink() or private_dir.stat().st_uid != os.getuid():
        raise PermissionError(f"{private_dir} is not this user's own directory")
    return private_dir / f"{run_dir.name}.series"


class Model:
    """Heliosphere 3D model output"""

    def __init__(self, dirname: Path, variable_mapping: dict):
        """
        Parameters
        ----------
        dirname : Path
            Path of the directory containing the model output
        variable_mapping : dict
            Mapping of variable names from app to model. This is
            helpful for keeping the app with a consistent naming
            while letting the models handle the variables themselves.
        """
        self.dir = dirname
        required_variables = {
            "velocity",
            "density",
            "pressure",
            "temperature",
            "b",
            "bx",
            "by",
            "bz",
            "dp",
        }
        for x in required_variables:
            if x not in variable_mapping:
                raise ValueError(
                    f"The required variable {x} was not found in the model's variable mapping."
                )
        self._variable_mapping = variable_mapping
        # Start off with a default Earth satellite filling the dictionary
        # keeping track of the satellites
        self.satellites = {"earth": ModelSatellite("earth", (-1, 0, 0))}

    def get_variable(self, name: str) -> str:
        """Get the variable name associated with this model

        Parameters
        ----------
        name : str
            The name of the variable within the App

        Returns
        -------
        str
            The name of the variable within the model
        """
        return self._variable_mapping[name]


class Enlil(Model):
    """
    The 3D Enlil model simulation results
    """

    def __init__(self, dirname: Path):
        variable_mapping = {
            "velocity": "Vr",
            "density": "Density",
            "pressure": "Pressure",
            "temperature": "T",
            "b": "Br",
            "bx": "Bx",
            "by": "By",
            "bz": "Bz",
            "dp": "DP",
        }
        super().__init__(dirname=dirname, variable_mapping=variable_mapping)

        self.data = pvs.NetCDFReader(
            registrationName="enlil-data", FileName=self._get_filenames(dirname)
        )
        self.data.Dimensions = ENLIL_DIMENSIONS
        self._load_satellites()

    @staticmethod
    def _get_filenames(dirname):
        """
        Get the file names to give the reader for a model run: always a
        single .series file. Once ParaView's file series reader has been
        given a .series file, it ignores a later plain list of file names
        (seen in ParaView 5.10.1, 6.1.1 and 6.2.0), so a change of run would
        show the previous run's data.

        There are two styles of run, old-style: one large file, new-style:
        individual file for each timestep.

        Parameters
        ----------
        dirname : Path
            Path of the directory containing the model output

        Returns
        -------
        list(str)
            The .series file name
        """
        legacy_filename = dirname / "pv-data-3d.nc"
        if legacy_filename.exists():
            # Old-style processing with a single giant file, which holds all
            # of the times itself
            series = private_series_path(dirname)
            write_series(series, [str(legacy_filename)])
            return [str(series)]

        # New processing with a single file for each timestep
        files = sorted(dirname.glob("pv-tim*.nc"))
        if not files:
            raise ValueError(f"No model output files in {dirname}")
        names = [f.name for f in files]
        # The run's own index, if it lists exactly these files
        series = dirname / SERIES_NAME
        if read_series(series) == names:
            return [str(series)]
        # Else one this server saved privately for exactly these files. A
        # private directory that is not this user's own, or one that cannot be
        # read (or vanished between checks), is skipped here, so the load can
        # still write the run directory's index below.
        try:
            private = private_series_path(dirname, create=False)
        except OSError:
            private = None
        if private and read_series(private, names_only=True) == [
            str(f) for f in files
        ]:
            return [str(private)]

        # No index of these files yet (or one of others): find each file's
        # time once and save them for later loads of this run
        times = file_times(files)
        if times is None:
            # The reader will have to ask each file for its times
            private = private_series_path(dirname)
            write_series(private, [str(f) for f in files])
            return [str(private)]
        try:
            write_series(series, names, times)
            return [str(series)]
        except OSError:
            # The run directory cannot be written to
            private = private_series_path(dirname)
            write_series(private, [str(f) for f in files], times)
            return [str(private)]

    def _load_satellites(self):
        sat_files = self.dir.glob("evo.*.json")
        for sat_file in sat_files:
            # strip the extra components from the name
            # Example: evo.sat_name.json
            name = sat_file.name[4:-5]
            self.satellites[name] = EnlilSatellite(name, sat_file)

    def change_run(self, dirname):
        """
        Change to a different model run.

        Parameters
        ----------
        dirname : Path
            Path of the directory containing the Enlil model output
        """
        # First, so a run with no files is refused before anything changes
        filenames = self._get_filenames(dirname)
        self.dir = dirname
        self.data.FileName = filenames
        # The reader reads a .series file again only when its name changes,
        # and this one may have been rewritten for files that arrived since
        self.data.SMProxy.UpdateProperty("FileName", 1)
        self.data.UpdatePipelineInformation()
        # Reload the satellite files for this run, starting again from the
        # default Earth as __init__ does, so a satellite this run has no file
        # for does not keep the previous run's track
        self.satellites = {"earth": ModelSatellite("earth", (-1, 0, 0))}
        self._load_satellites()


class Euhforia(Model):
    """
    The 3D EUHFORIA model simulation results
    """

    def __init__(self, dirname: Path):
        """
        Parameters
        ----------
        dirname : Path
            Path of the directory containing the EUHFORIA model output
        """
        variable_mapping = {
            "velocity": "vr",
            "density": "n-scaled",
            "pressure": "P",
            "temperature": "T",
            "b": "Br",
            "bx": "Bx",
            "by": "By",
            "bz": "Bz",
            "dp": "DP",
        }
        super().__init__(dirname=dirname, variable_mapping=variable_mapping)
        # Glob to list all files in the data directory, in name order as
        # change_run does
        fnames = sorted(str(fname) for fname in self.dir.glob("data_*.vts"))
        if not fnames:
            raise ValueError(f"No model output files in {self.dir}")
        # create a new 'XML Structured Grid Reader'
        self._input_data = pvs.XMLStructuredGridReader(
            registrationName="euhforia-data", FileName=fnames
        )
        self._input_data.CellArrayStatus = ["vr", "n", "P", "Br", "Bx", "By", "Bz"]
        self._input_data.TimeArray = "None"

        # Now scale all the arrays we need to work with
        self.data = pvs.CellDatatoPointData(
            registrationName="euhforia-pointdata", Input=self._input_data
        )
        self.data.ProcessAllArrays = 1
        self.data.PassCellData = 1

        # now calculate the density with rho * r**2
        self.data = pvs.Calculator(
            registrationName="euhforia-calculator", Input=self.data
        )
        self.data.AttributeType = "Point Data"
        self.data.ResultArrayName = "n-scaled"
        self.data.Function = "n * (coordsX^2 + coordsY^2 + coordsZ^2)"

        # Earth is at +1 X in Euhforia
        self.satellites["earth"].position = (1, 0, 0)

    def change_run(self, dirname):
        """
        Change to a different model run.

        Parameters
        ----------
        dirname : Path
            Path of the directory containing the EUHFORIA model output
        """
        # Refuse a run with no files before anything changes
        fnames = sorted(str(fname) for fname in dirname.glob("data_*.vts"))
        if not fnames:
            raise ValueError(f"No model output files in {dirname}")
        self.dir = dirname
        self._input_data.FileName = fnames
        # Read the new run's time steps now, as Enlil.change_run does, so the
        # caller's move to the first time step sees this run's times
        self._input_data.UpdatePipelineInformation()


class ModelSatellite:
    """
    A satellite or observation point in space.

    Parameters
    ----------
    name : str
        Name of the satellite
    """

    def __init__(self, name, position=(0, 0, 0)):
        self.name = name
        self.position = position

    def get_position(self, time):
        """The position of the satellite closest to the given time

        time : datetime-like
            Time of interest

        Returns
        -------
        The closest (X, Y, Z) position of the satellite to the requested time.
        """
        # Default just return the current position
        return self.position


class EnlilSatellite(ModelSatellite):
    def __init__(self, name, fpath):
        """
        Enlil time-series JSON of satellite data.

        This class will read in the given JSON file that contains
        the Enlil output for a specific satellite.

        Parameters
        ----------
        name : str
            Name of the satellite
        fpath : Path
            Path to the "evo" data file
        """
        super().__init__(name=name)
        with open(fpath) as f:
            self.json = json.loads(f.read())
            # all times are based off of unix epoch
            self.times = np.datetime64("1970-01-01") + np.array(
                self.json["coords"]["time"]["data"]
            ).astype("timedelta64[s]")

            # iterate over all data variables and set those as
            # attributes on the object
            for var in ["X", "Y", "Z"]:
                setattr(self, var, np.array(self.json["data_vars"][var]["data"]))

    def get_position(self, time):
        """
        Get the position at the given time.

        time : datetime-like
            Time of interest

        Returns
        -------
        Tuple (X, Y, Z) of the position of the satellite to the requested time.
        """
        loc = np.argmin(np.abs(np.datetime64(time) - self.times))
        return self.X[loc], self.Y[loc], self.Z[loc]
