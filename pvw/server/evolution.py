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
