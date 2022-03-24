import json
import os

import numpy as np


class Evolution:
    def __init__(self, fname):
        """
        Enlil time-series JSON dataset.

        This class will read in the given JSON file that contains
        the Enlil output for a specific satellite.
        """
        # fname will be something like: evo.earth.json
        # so, we will split on the periods and take the middle entry
        self.name = os.path.basename(fname).split('.')[1]
        with open(fname) as f:
            self.json = json.loads(f.read())
            # all times are based off of unix epoch
            self.times = (np.datetime64('1970-01-01') +
                          np.array(self.json['coords']['time']
                                   ['data']).astype('timedelta64[s]'))

            # iterate over all data variables and set those as
            # attributes on the object
            for var in ['X', 'Y', 'Z']:
                setattr(self, var,
                        np.array(self.json['data_vars'][var]['data']))

    def get_position(self, time):
        """
        Get the position at the given time.

        time : datetime-like
            Time of interest


        Returns
        -------
        The closest (X, Y, Z) position of the satellite to the requested time.
        """
        loc = np.argmin(np.abs(np.datetime64(time) - self.times))
        return self.X[loc], self.Y[loc], self.Z[loc]

    def get_data(self, variable):
        """
        Get the data within the variable of this satellite

        variable : str
            Variable of interest

        Returns
        -------
        List of data for this satellite
        """
        return self.json['data_vars'][variable]['data']

    def get_times(self):
        """
        Get the time series within the variable of this satellite

        Returns
        -------
        List of times for this satellite
        """
        return self.json['coords']['time']['data']

    def as_latis(self):
        """Create a Latis-style return for front-end use."""
        # List of lists
        # (ntimes, nvariables)
        ntimes = len(self.times)
        timestep_data = []
        for i in range(ntimes):
            curr_row = [self.get_times()[i]]
            for var in ["Density", "Vr", "Pressure", "T", "Bx", "By", "Bz"]:
                curr_row.append(self.get_data(var)[i][0])
            timestep_data.append(curr_row)

        json_out = {f"{self.name}": {
            "metadata": {
                "time": {
                    "units": "seconds since 1970-01-01",
                    "length": f"{ntimes}"
                },
                "density": {
                    "missing_value": "99999.99",
                    "description": "Density",
                    "units": "r<sup>2</sup>N/cm<sup>3</sup>"
                },
                "velocity": {
                    "missing_value": "99999.99",
                    "description": "Velocity",
                    "units": "km/s"
                },
                "pressure": {
                    "missing_value": "99999.99",
                    "description": "Ram pressure",
                    "units": ("r<sup>2</sup>N/cm<sup>3</sup> * "
                              "km<sup>2</sup>/s<sup>2</sup>")
                },
                "temperature": {
                    "missing_value": "99999.99",
                    "description": "Temperature",
                    "units": "K"
                },
                "bx": {
                    "missing_value": "99999.99",
                    "description": "BX",
                    "units": "nT"
                },
                "by": {
                    "missing_value": "99999.99",
                    "description": "BX",
                    "units": "nT"
                },
                "bz": {
                    "missing_value": "99999.99",
                    "description": "BX",
                    "units": "nT"
                }
            },
            "parameters": [
                "time",
                "density",
                "velocity",
                "pressure",
                "temperature",
                "bx",
                "by",
                "bz"
            ],
            "data": timestep_data
            }
        }
        return json.dumps(json_out)
