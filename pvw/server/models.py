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
