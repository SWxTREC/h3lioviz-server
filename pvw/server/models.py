from pathlib import Path

import paraview.simple as pvs


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
        self._data_dir = dirname
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
        fname = str(self._data_dir / "pv-data-3d.nc")
        self.data = pvs.NetCDFReader(registrationName="enlil-data", FileName=[fname])
        self.data.Dimensions = "(longitude, latitude, radius)"


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
        # Glob to list all files in the data directory
        fnames = [str(fname) for fname in self._data_dir.glob("data_*.vts")]
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
