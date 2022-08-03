from datetime import datetime, timedelta
import glob
import hashlib
import json
import os
import pathlib
import sys
import time
import urllib.request

import numpy as np
import xarray as xr

# Earth to Sun distance (m)
AU = 1.496e11
# Mass of hydrogen (kg)
m_hydrogen = 1.6735575e-27
# m3 to cm3
m3_to_cm3 = (1.0 / 100) ** 3
# kg/m3 -> N/cm3
density_conversion = 1.0 / m_hydrogen * m3_to_cm3
# m/s -> km/s
velocity_conversion = 1.0 / 1000.0


class NumpyEncoder(json.JSONEncoder):
    """Convert numpy int/float/array to JSON serializable objects."""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def spherical_to_cartesian(r, theta, phi, v1, v2, v3):
    """
    Transform vector field components from spherical to cartesian.

    The original coordinates are HEEQ+180, keep that system for now
    and put the X/Y/Z in that frame.
    """
    # v1, v2, v3 are the vector components in the three directions
    vx = (
        v1 * np.sin(theta) * np.cos(phi)
        + v2 * np.cos(theta) * np.cos(phi)
        - v3 * np.sin(phi)
    )
    vy = (
        v1 * np.sin(theta) * np.sin(phi)
        + v2 * np.cos(theta) * np.sin(phi)
        + v3 * np.cos(phi)
    )
    vz = v1 * np.cos(theta) - v2 * np.sin(theta)

    return vx, vy, vz


def process_tim(ds):
    """
    Function to process each tim file.

    This does coordinate transformations and renaming.
    """
    ds = ds[{"nblk": 0}]
    ds["n1"] = ds["X1"]
    ds["n2"] = ds["X2"]
    ds["n3"] = ds["X3"]

    try:
        t0 = datetime.strptime(ds.attrs["rundate_cal"], "%Y-%m-%dT%H")
    except ValueError:
        # May not be an hour given?
        t0 = datetime.strptime(ds.attrs["rundate_cal"], "%Y-%m-%d")
    t = t0 + timedelta(seconds=ds["TIME"].item())

    # Change from Tesla to nT
    for x in ["B1", "B2", "B3"]:
        # Multiply by the polarity tracer to get the proper direction (+/-)
        # The polarity is just a "tracer" variable with no physical meaning
        # it is only used to get the sign of the magnetic field right
        ds[x] *= 1e9 * np.sign(ds["BP"])

    # Br, Blat, Blon
    ds["Bx"], ds["By"], ds["Bz"] = spherical_to_cartesian(
        ds["n1"], ds["n2"], ds["n3"], ds["B1"], ds["B2"], ds["B3"]
    )
    ds["Br"] = ds["B1"]

    # ds['Vx'], ds['Vy'], ds['Vz'] =
    # spherical_to_cartesian(ds['n1'], ds['n2'], ds['n3'],
    # ds['V1'], ds['V2'], ds['V3'])

    name_map = {
        "n1": "radius",
        "n2": "latitude",
        "n3": "longitude",
        "D": "Density",
        "TIME": "time",
    }

    ds = ds.rename(name_map)
    ds["radius"] = ds["radius"] / AU
    ds["latitude"] = np.rad2deg(np.pi / 2 - ds["latitude"])
    # subtract pi to point towards Earth
    ds["longitude"] = np.rad2deg(ds["longitude"])

    ds["Vr"] = ds["V1"] * velocity_conversion
    ds["Density"] = (
        ds["Density"]
        * density_conversion
        / (1 - ds.attrs["xalpha"])
        * ds["radius"] ** 2
    )
    ds["DP"] *= ds["radius"] ** 2
    # Ram pressure (rho * v**2)
    ds["Pressure"] = ds["Density"] * ds["Vr"] ** 2

    ds["radius"].attrs.update({"units": "AU"})
    ds["latitude"].attrs.update({"units": "degrees_north"})
    ds["longitude"].attrs.update({"units": "degrees_east"})
    ds["time"] = t
    # TODO: Do we want to be able to display a string in Paraview?
    #       I'm not sure if we should convert this to string here,
    #       or use the integers and index within Paraview.
    # ds['time_string'] = datetime.strftime(t, '%Y-%m-%dT%H:00')
    # ds['date'] = xr.DataArray([t.year, t.month, t.hour, t.minute])
    # ds['time'] = ds['time'].dt.round('min') # Strip the seconds off
    ds = ds.assign_coords({"time": ds["time"]}).expand_dims("time")

    # Remove unnecessary variables to make file sizes smaller
    ds = ds.drop_vars(
        ["DT", "NSTEP", "BP", "X1", "X2", "X3", "B1", "B2", "B3", "V1", "V2", "V3"]
    )

    return ds


def process_evo(ds):
    """
    Function to process each evo file.

    This does coordinate transformations and renaming.
    """
    try:
        t0 = datetime.strptime(ds.attrs["rundate_cal"], "%Y-%m-%dT%H")
    except ValueError:
        # May not be an hour given?
        t0 = datetime.strptime(ds.attrs["rundate_cal"], "%Y-%m-%d")
    t = np.datetime64(t0) + np.timedelta64(1, "s") * ds["TIME"]

    # Change from Tesla to nT
    for x in ["B1", "B2", "B3"]:
        # Multiply by the polarity tracer to get the proper direction (+/-)
        # The polarity is just a "tracer" variable with no physical meaning
        # it is only used to get the sign of the magnetic field right
        ds[x] *= 1e9 * np.sign(ds["BP"])

    # Br, Blat, Blon
    ds["Bx"], ds["By"], ds["Bz"] = spherical_to_cartesian(
        ds["X1"], ds["X2"], ds["X3"], ds["B1"], ds["B2"], ds["B3"]
    )

    # ds['Vx'], ds['Vy'], ds['Vz'] =
    # spherical_to_cartesian(ds['n1'], ds['n2'], ds['n3'],
    # ds['V1'], ds['V2'], ds['V3'])

    name_map = {"D": "Density"}

    ds = ds.rename(name_map)
    rad = ds["X1"] / AU
    lat = np.pi / 2 - ds["X2"]
    # subtract pi to point towards Earth
    lon = ds["X3"]
    ds["X"] = rad * np.cos(lon) * np.cos(lat)
    ds["Y"] = rad * np.sin(lon) * np.cos(lat)
    ds["Z"] = rad * np.sin(lat)

    ds["Vr"] = ds["V1"] * velocity_conversion
    ds["Density"] = (
        ds["Density"] * density_conversion / (1 - ds.attrs["xalpha"]) * rad ** 2
    )
    # Ram pressure (rho * v**2)
    ds["Pressure"] = ds["Density"] * ds["Vr"] ** 2

    ds["X"].attrs.update({"units": "AU"})
    ds["Y"].attrs.update({"units": "AU"})
    ds["Z"].attrs.update({"units": "AU"})
    ds = ds.rename_dims({"nevo": "time"})
    ds = ds.assign_coords(time=t.data)

    # Remove unnecessary variables to make file sizes smaller
    ds = ds.drop_vars(
        [
            "DT",
            "NSTEP",
            "BP",
            "TIME",
            "X1",
            "X2",
            "X3",
            "B1",
            "B2",
            "B3",
            "V1",
            "V2",
            "V3",
        ]
    )

    return ds


def process_directory(path):
    print("Beginning processing, may take several minutes")
    t0 = time.time()

    # ---------
    # TIM file processing
    # ---------
    fnames = sorted(glob.glob(path + "/tim.*.nc"))
    if len(fnames) == 0:
        raise ValueError("No files found to process in the current directory.")

    # Load into a multi-file dataset to be able to concatenate along Time
    ds = xr.open_mfdataset(
        fnames, combine="by_coords", preprocess=process_tim, engine="netcdf4"
    )
    print(f"Dataset loaded: {time.time()-t0} s")
    # New path based on the metadata
    s = json.dumps(ds.attrs, cls=NumpyEncoder)

    run_id = hashlib.sha256(s.encode("utf-8")).hexdigest()[:8]
    newpath = path + "/pv-ready-data-" + run_id
    if not os.path.exists(newpath):
        os.mkdir(newpath)

    with open(f"{newpath}/metadata.json", "w") as f:
        # s is currently just the string version of the json dict
        # because we can't hash a dictionary. So, lets unpack s
        # back to a dict and add the run_id to it for later reference.
        d = json.loads(s)
        d["run_id"] = run_id
        f.write(json.dumps(d))

    ds.to_netcdf(
        f"{newpath}/pv-data-3d.nc",
        engine="scipy",
        encoding={"time": {"units": "seconds since 1970-01-01"}},
    )
    # _, datasets = zip(*ds.groupby('time'))
    # paths = [f"{newpath}/tim.{i:04d}.nc" for i in range(len(fnames))]
    # xr.save_mfdataset(datasets, paths, engine='scipy')
    print(f"Dataset saved: {time.time()-t0} s")

    # Load and process the evolution (evo) file
    fnames = sorted(glob.glob(path + "/evo.*.nc"))
    if len(fnames) == 0:
        raise ValueError(
            "No evolution files found to process in the " "current directory."
        )

    datasets = [
        process_evo(xr.open_dataset(fname, engine="netcdf4")) for fname in fnames
    ]

    for ds_evo, fname in zip(datasets, fnames):
        # Save a single evolution file to NetCDF
        newfile = f"{newpath}/{os.path.basename(fname)}"
        ds_evo.to_netcdf(
            newfile, encoding={"time": {"units": "seconds since 1970-01-01"}}
        )
        # Reload that file and save it to json
        # Loading with decode_times=False makes it so the datetimes
        # get encoded as ints within the JSON
        newds = xr.load_dataset(newfile, decode_times=False)

        # Drop the ".nc" extension and replace with json
        with open(newfile.replace(".nc", ".json"), "w") as f:
            # Put it through dump/load/dump to get floating point truncation
            f.write(
                json.dumps(
                    json.loads(
                        json.dumps(newds.to_dict()),
                        parse_float=lambda x: round(float(x), 3),
                    )
                )
            )

    print(f"Evo datasets saved: {time.time()-t0} s")

    # Downloading images now, we want to download for every day in the dataset
    # Convert numpy datetime64 (strip nanoseconds component),
    # then to a Python datetime object
    start_date = ds["time"].data[0].astype("datetime64[s]").astype(object)
    end_date = ds["time"].data[-1].astype("datetime64[s]").astype(object)
    dt = timedelta(days=1)
    curr_date = start_date
    while curr_date <= end_date:
        download_hmi(curr_date, outdir=newpath + "/solar_images")
        curr_date += dt

    print(f"Images saved: {time.time()-t0} s")


def download_hmi(date: datetime, outdir=None, resolution="1k"):
    """
    Downloads hourly HMI files for the date selected if not already
    available locally.
    """
    if outdir is None:
        outdir = pathlib.Path.cwd() / "solar_images"
    else:
        outdir = pathlib.Path(outdir)

    if not outdir.exists():
        # Make our directory structure if it doesn't exist yet
        print(f"Making directory: {outdir}")
        outdir.mkdir(parents=True, exist_ok=True)

    base_url = "http://jsoc.stanford.edu/data/hmi/images/"
    # Files are stored in YYYY/MM/DD directories, with
    # specific times in the filenames after that.
    # This is the alternative request to search for nearest that gets
    # forwarded to the below url
    # ("http://jsoc.stanford.edu/cgi-bin/hmiimage.pl"
    #        "?Year=2017&Month=09&Day=06&Hour=00&Minute=00"
    #        "&Kind=_M_color_&resolution=4k")
    # url = ("http://jsoc.stanford.edu/data/hmi/images/"
    #        "2017/09/06/20170906_000000_M_color_4k.jpg")
    dt = timedelta(hours=1)
    for timestep in range(24):
        t = date + timestep * dt
        url_dir = t.strftime("%Y/%m/%d/")
        fname = t.strftime(f"%Y%m%d_%H0000_M_color_{resolution}.jpg")
        url = base_url + url_dir + fname
        download_path = outdir / fname
        if download_path.exists():
            # We have already downloaded this file, so no need to repeat
            print(f"Already downloaded, skipping: {fname}")
            continue

        # Make the download request
        try:
            req = urllib.request.urlopen(url)
        except urllib.error.HTTPError:
            # Just ignore the URL not found errors for now
            print(f"Could not download {url}")
            continue

        # Write out the response to our local file
        with open(download_path, "wb") as f:
            f.write(req.read())
        print(f"Downloaded {url}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise ValueError(
            "This program takes one argument, the path to the " "directory to convert"
        )
    process_directory(sys.argv[1])
