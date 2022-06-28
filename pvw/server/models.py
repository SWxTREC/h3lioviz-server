from pathlib import Path

import paraview.simple as pvs


class Enlil:
    """
    The 3D Enlil model simulation results
    """

    def __init__(self, dirname: Path):
        """
        Parameters
        ----------
        dirname : Path
            Path of the directory containing the Enlil model output
        """
        self._data_dir = dirname
        fname = str(self._data_dir / "pv-data-3d.nc")
        self.data = pvs.NetCDFReader(registrationName="enlil-data", FileName=[fname])
        self.data.Dimensions = "(longitude, latitude, radius)"

        # Use properties for the variable mapping
        var_map = {
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
        for key, val in var_map.items():
            setattr(self, key, val)


class Euhforia:
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
        self._data_dir = dirname
        fnames = [str(fname) for fname in dirname.glob("data_*.vts")]
        # create a new 'XML Structured Grid Reader'
        self.data = pvs.XMLStructuredGridReader(
            registrationName="euhforia-data", FileName=fnames
        )
        self.data.CellArrayStatus = ["vr", "n", "P", "Br", "Bx", "By", "Bz"]
        self.data.TimeArray = "None"

        # Now scale all the arrays we need to work with
        self.data = pvs.CellDatatoPointData(
            registrationName="euhforia-pointdata", Input=self.data
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

        # Use properties for the variable mapping
        var_map = {
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
        for key, val in var_map.items():
            setattr(self, key, val)
