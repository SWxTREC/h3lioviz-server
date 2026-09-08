import json
import pathlib
import re
from datetime import datetime, timedelta, timezone

import numpy as np
import requests
import xarray as xr

# Converts names of objects to the codes needed to request position data from NASA HELIOWeb
HELIOWEB_OBJECT_CODES = {
    "earth": "04",
    "mars": "42",
    "venus": "15",
    "mercury": "29",
    "solar_orbiter": "46",
}
J2000_JD = 2451545.0
# Number of days in a julian century
JULIAN_CENTURY_DAYS = 36525.0


def process_helioweb(
    satellite_list: list[str],
    start_date: datetime,
    end_date: datetime,
    newpath: pathlib.Path,
):
    for satellite in satellite_list:
        # Download NASA helioweb data for each satellite
        hgi_data = download_helioweb_hgi(satellite, start_date, end_date)

        # Convert from NASA helioweb's coordinates, HGI,
        # to the system used by ENLIL, HNM (HEEQ+180).
        hnm_data = hgi_to_hnm_from_date(
            hgi_data["year"],
            hgi_data["day"],
            hgi_data["hour"],
            hgi_data["radius_au"],
            hgi_data["hgi_lat"],
            hgi_data["hgi_lon"],
        )

        # Convert from ENLIL's HNM to what paraview expects.
        # This means converting from spherical to cartesian
        # as well as setting up time as seconds since 1970-01-01
        pv_data = hnm_to_paraview(
            hnm_data["year"],
            hnm_data["day"],
            hnm_data["hour"],
            hnm_data["radius_au"],
            hnm_data["colattitude_rad"],
            hnm_data["longitude_rad"],
        )

        # Build the processed xr.Dataset that will be converted
        # to both a netCDF file # and json file for use by paraview.
        # All of the variables besides X, Y, Z, and time are set to NaN
        ds = build_pv_evo_dataset(
            pv_data["X"], pv_data["Y"], pv_data["Z"], pv_data["time"]
        )

        # TODO we could probably have a shared "save_evo()" function for both helioweb
        # and normal ENLIL evo.
        # Copied from the evo logic in process_directory()
        newfile = newpath / f"evo.{satellite}.nc"
        ds.to_netcdf(newfile)

        newds = xr.load_dataset(newfile, decode_times=False)

        json_file = pathlib.Path(str(newfile).replace(".nc", ".json"))
        with open(json_file, "w") as f:
            f.write(
                json.dumps(
                    json.loads(
                        json.dumps(newds.to_dict()),
                        parse_float=lambda x: round(float(x), 3),
                    )
                )
            )


def download_helioweb_hgi(object_name: str, start_date: datetime, end_date: datetime):
    """
    Download hourly HGI coordinates from NASA HELIOWeb for a given object
    Web GUI: https://omniweb.gsfc.nasa.gov/coho/helios/heli.html
    """
    query = {
        "activity": "retrieve",
        "object": HELIOWEB_OBJECT_CODES[object_name],
        # This cooresponds to the HGI coordinate system
        "coordinate": "3",
        "resolution": "Hourly",
        # Number of decimal places for radius (needs verification maybe it applies to lat/lon as well
        "precn": "4",
        "start_year": start_date.year,
        "start_day": start_date.timetuple().tm_yday,
        "stop_year": end_date.year,
        "stop_day": end_date.timetuple().tm_yday,
        # Determines which equinox epoch to use. This corresponds to which instantaneous "First Point of Aries" (FPA) direction is used to specify longitude
        # 2 corresponds to Mean-of-Date which is a smoothed FPA at 00.00 on date of interest
        "equinox": "2",
    }

    # Get the data from HELIOWeb
    response = requests.post(
        "https://omniweb.gsfc.nasa.gov/cgi/models/helios2_h.cgi", data=query, timeout=30
    )
    # Raises an error on any non 2XX status code
    response.raise_for_status()

    # Extract the table from the html
    text = response.text
    match = re.search(r"<pre>(.*?)</pre>", text, flags=re.DOTALL)
    if not match:
        raise ValueError(
            f"The response from HELIOWeb for {object_name} did not match as expected."
        )
    table_text = match.group(1)

    # Convert the table text into a list of tuples which can then be converted into a np array
    rows = []
    for line in table_text.splitlines():
        # Each line should just be whitespace separated text
        cols = line.split()
        # Skip the table header
        if cols[0] == "YEAR":
            continue

        if len(cols) != 6 or not cols[0].isdigit():
            print(line)
            print(
                f"Warning: malformed line in downloaded HELIOWeb data for {object_name}. See above output for details. Continuing..."
            )
            continue

        rows.append(tuple(cols))

    output = np.array(
        rows,
        dtype=[
            ("year", "i4"),
            ("day", "i4"),
            ("hour", "i4"),
            ("radius_au", "f8"),
            ("hgi_lat", "f8"),
            ("hgi_lon", "f8"),
        ],
    )

    return output


def hgi_to_hnm_from_date(year, day_of_year, hour, r_hgi, lat_hgi, lon_hgi):
    """Convert HGI Spherical coordinates to HNM spherical"""
    jd = julian_date(year, day_of_year, hour)

    lambda_angle = lambda_from_julian_date(jd)
    omega_angle = omega_from_julian_date(jd)

    rotation_angle = lambda_angle - omega_angle

    # shift by rotation angle for HGI -> HEEQ
    lon_heeq = np.deg2rad(lon_hgi) - rotation_angle
    # 180 degree shift for HEEQ -> HNM
    lon_hnm = lon_heeq + np.pi
    lon_hnm = normalize_radians(lon_hnm)

    # Latitude defines 90 degrees as the +Z axis and -90 as the -Z axis
    # Co-latitude defines 0 degrees as the +Z axis and 180 as the -Z axis
    # The coordinate system for ENLIL uses colatitude so we convert
    colat_hnm = normalize_radians(-(np.deg2rad(lat_hgi) - np.pi / 2))

    # There are a few differences between the outputted numpy array and the ENLIL .nc files
    # The ENLIL .nc files have time measured in seconds relative to rundate_cal
    # The ENLIL .nc files use radius in meters
    data = np.array(
        list(zip(year, day_of_year, hour, r_hgi, colat_hnm, lon_hnm)),
        dtype=[
            (
                "year",
                "i4",
            ),
            ("day", "i4"),
            ("hour", "i4"),
            ("radius_au", "f8"),
            ("colattitude_rad", "f8"),
            ("longitude_rad", "f8"),
        ],
    )
    return data


def t0_from_julian_date(jd):
    """t0 is Julian centuries from J2000.0."""
    return (np.asarray(jd) - J2000_JD) / JULIAN_CENTURY_DAYS


def lambda_from_julian_date(jd):
    """
    Approximate heliocentric ecliptic longitude of Earth at J2000.

    This uses the standard low-precision solar ephemeris: calculate the Sun's
    true geocentric ecliptic longitude, then add 180 degrees for Earth's
    heliocentric longitude. The result is relative to the ecliptic of date.
    """
    # Number of centuries since J2000
    t0 = t0_from_julian_date(jd)
    # Mean geometric geocentric ecliptic longitude of the sun.
    # Geometric, as opposed to apparent, refers to the Sun's physical location and not the location in the sky with the effects of our atmosphere
    # Source of equation: 'Geom Mean Long Sun' https://gml.noaa.gov/grad/solcalc/NOAA_Solar_Calculations_day.xls
    mean_longitude = 280.46646 + 36000.76983 * t0 + 0.0003032 * t0**2
    # Mean anomaly is the fraction of Earth's idealized circular orbit that has elapsed from the perihelion
    # expressed as an angle for a given time. The first term, 357.52911, is the Earth's
    # mean anomaly at J2000. Earth advances 359.990503° in mean anomaly during one Julian year (365.25 days),
    # measured relative to the moving perihelion direction and is called the mean anomaly rate which
    # is an angular velocity. Although this mean anomaly rate is assumed to be fixed for a given year,
    # it drifts with time due to a variety of effects.
    # Source of equation: 'Geom Mean Anom Sun' https://gml.noaa.gov/grad/solcalc/NOAA_Solar_Calculations_day.xls
    mean_anomaly = np.deg2rad(357.52911 + 35999.05029 * t0 - 0.0001537 * t0**2)
    # Equation of the center is the correction needed from mean anomaly to true anomaly.
    # In other words it is the difference between the true anomaly and the mean anomaly.
    # Because the earth is a low-ecentricity orbit that varies over time, we can model
    # it with a time varrying sine series.
    # We should expect the equation_of_center to be 0 at the perihelion and apehelion
    # As the mean anomaly and true anomaly orbits intersect at these points. This correction for the Earth
    # should max out at around 2 degrees for a given time.
    # Source of equation: 'Sun Eq of Ctr' https://gml.noaa.gov/grad/solcalc/NOAA_Solar_Calculations_day.xls
    equation_of_center = (
        (1.914602 - 0.004817 * t0 - 0.000014 * t0**2) * np.sin(mean_anomaly)
        + (0.019993 - 0.000101 * t0) * np.sin(2.0 * mean_anomaly)
        + 0.000289 * np.sin(3.0 * mean_anomaly)
    )
    # Get the true heliocentric ecliptic longitude of the earth by
    # appplying the correction, equation of the center, to the mean geocentric ecliptic longitude of the sun
    # and apply a 180 degree rotation to get the mean heliocentric ecliptic longitude of the earth.
    # This is just a rotation about the Z-axis, solar rotation axis to face the Earth rather than the Sun.
    # We just need to add the equation of the center (difference between the mean anomaly and true anomaly):
    # equation_of_center = true_anomaly - mean_anomaly
    # true_anomaly = equation_of_center + mean_anomaly
    # Lmean = Lperihelion + mean_anomaly
    # Ltrue = Lperihilion + true_anomaly
    # Ltrue = Lperihilion + mean_anomaly + equation_of_center
    # Ltrue = Lmean + equation_of_center
    # Source of equation: 'Sun True Long' https://gml.noaa.gov/grad/solcalc/NOAA_Solar_Calculations_day.xls
    lambda_angle = mean_longitude + equation_of_center + 180.0
    return np.deg2rad(lambda_angle)


def julian_date(year, day_of_year, hour=0.0):
    """
    Convert year, day-of-year, and decimal hour to Julian Date.

    Inputs may be scalars or array-like values. Day-of-year is one-indexed,
    matching the values in the HGI .lst files.
    """
    year_arr, day_arr, hour_arr = np.broadcast_arrays(year, day_of_year, hour)

    def scalar_julian_date(y, d, h):
        start = datetime(int(y), 1, 1, tzinfo=timezone.utc)
        dt = start + timedelta(days=float(d) - 1.0, hours=float(h))
        unix_days = dt.timestamp() / 86400.0
        return 2440587.5 + unix_days

    if year_arr.shape == ():
        return scalar_julian_date(year_arr.item(), day_arr.item(), hour_arr.item())

    vectorized = np.vectorize(scalar_julian_date, otypes=[float])
    return vectorized(year_arr, day_arr, hour_arr)


def omega_from_julian_date(jd):
    """
    Longitude of the ascending node of the solar equator.
    Curl your right hand in the direction of earth's orbit (counter clockwise)
    and the ascending node is where the ecliptic

    The traditional definition is relative to the ecliptic plane of date:
    Omega = 75.76 deg + 1.397 deg * T0.
    """
    t0 = t0_from_julian_date(jd)
    omega = 75.76 + 1.397 * np.asarray(t0)

    return np.deg2rad(omega)


def hnm_to_paraview(year, day_of_year, hour, r, colat, lon):
    # Do a standard spherical to cartesian conversion
    x = r * np.sin(colat) * np.cos(lon)
    y = r * np.sin(colat) * np.sin(lon)
    z = r * np.cos(colat)

    # Convert year, day_of_year, hour to seconds since 1970-01-01
    times = np.array(
        [
            (
                datetime(int(year), 1, 1, tzinfo=timezone.utc)
                + timedelta(days=int(day) - 1, hours=int(hour))
            ).timestamp()
            for year, day, hour in zip(year, day_of_year, hour)
        ]
    )

    data = np.asarray(
        list(zip(x, y, z, times)),
        dtype=[("X", "f8"), ("Y", "f8"), ("Z", "f8"), ("time", "i4")],
    )

    return data


def build_pv_evo_dataset(x, y, z, t):

    shape = len(t)
    nan = np.full(shape, np.nan, dtype=np.float64)
    nanf = np.full(shape, np.nan, dtype=np.float32)

    # The decision to use nan/nanf was determined by
    # looking at the datatypes in the processed evo.earth.nc
    # using `ncdump -h <path-to-processed-evo-earth-nc>`
    ds = xr.Dataset(
        data_vars={
            "X": xr.Variable("time", x, attrs={"units": "AU"}),
            "Y": xr.Variable("time", y, attrs={"units": "AU"}),
            "Z": xr.Variable("time", z, attrs={"units": "AU"}),
            "Density": xr.Variable("time", nan),
            "T": xr.Variable(
                "time", nanf, attrs={"units": "K", "long_name": "Temperature"}
            ),
            # The current pv-evo datasets don't do a great job of defining these fields so not adding additional
            # attrs for now
            # TODO explore these more and provide proper units + descriptions
            "DP": xr.Variable("time", nanf),
            "Bx": xr.Variable("time", nanf),
            "By": xr.Variable("time", nanf),
            "Bz": xr.Variable("time", nanf),
            "Br": xr.Variable("time", nanf),
            "Vr": xr.Variable("time", nanf),
            "Pressure": xr.Variable("time", nan),
        },
        coords={
            "time": xr.Variable(
                "time",
                t,
                attrs={
                    "units": "seconds since 1970-01-01",
                    "calendar": "proleptic_gregorian",
                },
            )
        },
    )

    return ds


def normalize_radians(angle):
    """Normalize angle(s) to [0, 2*pi) radians."""
    return np.mod(angle, 2.0 * np.pi)
