import io
import json
import os

import boto3

from helpers import _get_env_var


class Evolution:
    def __init__(self, name, json_data):
        """
        Enlil time-series JSON dataset.

        This class will read in the given JSON file that contains
        the Enlil output for a specific satellite.
        """
        self.name = name
        self.json = json_data

    def get_data(self, variable):
        """
        Get the data within the variable of this satellite

        variable : str
            Variable of interest

        Returns
        -------
        List of data for this satellite
        """
        return self.json["data_vars"][variable]["data"]

    def get_times(self):
        """
        Get the time series within the variable of this satellite

        Returns
        -------
        List of times for this satellite
        """
        return self.json["coords"]["time"]["data"]

    def as_latis(self):
        """Create a Latis-style return for front-end use."""
        # List of lists
        # (ntimes, nvariables)
        times = self.get_times()
        ntimes = len(times)
        timestep_data = []
        for i in range(ntimes):
            # Time is in units of seconds from epoch, but Scicharts needs
            # the units to be in milliseconds
            curr_row = [int(times[i]) * 1000]
            for var in ["Density", "Vr", "Pressure", "T", "Bx", "By", "Bz"]:
                curr_data = self.get_data(var)[i]
                # The data can be multi-block, which would mean this is a list,
                # and we want the first value out of the list from block-0, which
                # is the base run. The other blocks are for separate CMEs
                if isinstance(curr_data, list):
                    curr_data = curr_data[0]
                curr_row.append(curr_data)
            timestep_data.append(curr_row)

        json_out = {
            f"{self.name}": {
                "metadata": {
                    "time": {
                        "units": "milliseconds since 1970-01-01",
                        "length": f"{ntimes}",
                    },
                    "density": {
                        "missing_value": "99999.99",
                        "description": "Density",
                        "units": "r<sup>2</sup>N/cm<sup>3</sup>",
                    },
                    "velocity": {
                        "missing_value": "99999.99",
                        "description": "Velocity",
                        "units": "km/s",
                    },
                    "pressure": {
                        "missing_value": "99999.99",
                        "description": "Ram pressure",
                        "units": (
                            "r<sup>2</sup>N/cm<sup>3</sup> * "
                            "km<sup>2</sup>/s<sup>2</sup>"
                        ),
                    },
                    "temperature": {
                        "missing_value": "99999.99",
                        "description": "Temperature",
                        "units": "K",
                    },
                    "bx": {
                        "missing_value": "99999.99",
                        "description": "BX",
                        "units": "nT",
                    },
                    "by": {
                        "missing_value": "99999.99",
                        "description": "BY",
                        "units": "nT",
                    },
                    "bz": {
                        "missing_value": "99999.99",
                        "description": "BZ",
                        "units": "nT",
                    },
                },
                "parameters": [
                    "time",
                    "density",
                    "velocity",
                    "pressure",
                    "temperature",
                    "bx",
                    "by",
                    "bz",
                ],
                "data": timestep_data,
            }
        }
        return json_out


def get_timeseries(run_id, satellite):
    bucket = _get_env_var("S3_BUCKET_NAME")

    s3_client = boto3.client("s3")

    # In case the frontend passes .jsond in, turn it into a trailing .json
    satellite = satellite.rstrip("d")
    key = f"data/h3lioviz/pv-ready-data-{run_id}/evo.{satellite}"
    # Create an in-memory buffer to read our returns into
    buffer = io.BytesIO()
    s3_client.download_fileobj(Bucket=bucket, Key=key, Fileobj=buffer)
    # Load the json metadata and append it to our content list
    evo = Evolution(satellite.replace(".json", ""), json.loads(buffer.getvalue()))

    response = {
        "statusCode": 200,
        "body": json.dumps(evo.as_latis()),
        "headers": {
            "Content-Type": "application/json",
        },
    }

    return response
