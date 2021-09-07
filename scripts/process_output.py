from datetime import datetime, timedelta
import glob
import json
import os
import sys
import time

import numpy as np
import xarray as xr

# Earth to Sun distance (m)
AU = 1.496e+11
# Mass of hydrogen (kg)
m_hydrogen = 1.6735575e-27
# m3 to cm3
m3_to_cm3 = (1./100)**3
# kg/m3 -> N/cm3
density_conversion = 1./m_hydrogen * m3_to_cm3
# m/s -> km/s
velocity_conversion = 1./1000.


def spherical_to_cartesian(r, theta, phi, v1, v2, v3):
    """
    Transform vector field components from spherical to cartesian.

    The original coordinates are HEEQ+180, keep that system for now
    and put the X/Y/Z in that frame.
    """
    # v1, v2, v3 are the vector components in the three directions
    vx = v1*np.sin(theta)*np.cos(phi) + v2*np.cos(theta) * \
        np.cos(phi) - v3*np.sin(phi)
    vy = v1*np.sin(theta)*np.sin(phi) + v2*np.cos(theta) * \
        np.sin(phi) + v3*np.cos(phi)
    vz = v1*np.cos(theta) - v2*np.sin(theta)

    return vx, vy, vz


def process_tim(ds):
    """
    Function to process each tim file.

    This does coordinate transformations and renaming.
    """
    ds = ds[{'nblk': 0}]
    ds['n1'] = ds['X1']
    ds['n2'] = ds['X2']
    ds['n3'] = ds['X3']

    t0 = datetime.strptime(ds.attrs['rundate_cal'], '%Y-%m-%dT%H')
    t = t0 + timedelta(seconds=ds['TIME'].item())

    # Change from Tesla to nT
    for x in ['B1', 'B2', 'B3']:
        # Multiply by the polarity tracer to get the proper direction (+/-)
        # The polarity is just a "tracer" variable with no physical meaning
        # it is only used to get the sign of the magnetic field right
        ds[x] *= 1e9 * np.sign(ds['BP'])

    # Br, Blat, Blon
    ds['Bx'], ds['By'], ds['Bz'] = spherical_to_cartesian(
        ds['n1'], ds['n2'], ds['n3'], ds['B1'], ds['B2'], ds['B3'])

    # ds['Vx'], ds['Vy'], ds['Vz'] =
    # spherical_to_cartesian(ds['n1'], ds['n2'], ds['n3'],
    # ds['V1'], ds['V2'], ds['V3'])

    name_map = {'n1': 'radius', 'n2': 'latitude', 'n3': 'longitude',
                'D': 'Density',
                'TIME': 'time'}

    ds = ds.rename(name_map)
    ds['radius'] = ds['radius']/AU
    ds['latitude'] = np.rad2deg(np.pi/2 - ds['latitude'])
    # subtract pi to point towards Earth
    ds['longitude'] = np.rad2deg(ds['longitude'])

    ds['Vr'] = ds['V1'] * velocity_conversion
    ds['Density'] = ds['Density'] * density_conversion / \
        (1 - ds.attrs['xalpha']) * ds['radius']**2
    # Ram pressure (rho * v**2)
    ds['Pressure'] = ds['Density'] * ds['Vr']**2

    ds['radius'].attrs.update({'units': 'AU'})
    ds['latitude'].attrs.update({'units': 'degrees_north'})
    ds['longitude'].attrs.update({'units': 'degrees_east'})
    ds['time'] = t
    # TODO: Do we want to be able to display a string in Paraview?
    #       I'm not sure if we should convert this to string here,
    #       or use the integers and index within Paraview.
    # ds['time_string'] = datetime.strftime(t, '%Y-%m-%dT%H:00')
    # ds['date'] = xr.DataArray([t.year, t.month, t.hour, t.minute])
    # ds['time'] = ds['time'].dt.round('min') # Strip the seconds off
    ds = ds.assign_coords({'time': ds['time']}).expand_dims('time')

    # Remove unnecessary variables to make file sizes smaller
    ds = ds.drop_vars(['DT', 'NSTEP', 'BP',
                       'X1', 'X2', 'X3',
                       'B1', 'B2', 'B3',
                       'V1', 'V2', 'V3'])

    return ds


def process_evo(ds):
    """
    Function to process each evo file.

    This does coordinate transformations and renaming.
    """
    t0 = datetime.strptime(ds.attrs['rundate_cal'], '%Y-%m-%dT%H')
    t = np.datetime64(t0) + np.timedelta64(1, 's') * ds['TIME']

    # Change from Tesla to nT
    for x in ['B1', 'B2', 'B3']:
        # Multiply by the polarity tracer to get the proper direction (+/-)
        # The polarity is just a "tracer" variable with no physical meaning
        # it is only used to get the sign of the magnetic field right
        ds[x] *= 1e9 * np.sign(ds['BP'])

    # Br, Blat, Blon
    ds['Bx'], ds['By'], ds['Bz'] = spherical_to_cartesian(
        ds['X1'], ds['X2'], ds['X3'], ds['B1'], ds['B2'], ds['B3'])

    # ds['Vx'], ds['Vy'], ds['Vz'] =
    # spherical_to_cartesian(ds['n1'], ds['n2'], ds['n3'],
    # ds['V1'], ds['V2'], ds['V3'])

    name_map = {'D': 'Density'}

    ds = ds.rename(name_map)
    rad = ds['X1']/AU
    lat = np.pi/2 - ds['X2']
    # subtract pi to point towards Earth
    lon = ds['X3']
    ds['X'] = rad * np.cos(lon) * np.cos(lat)
    ds['Y'] = rad * np.sin(lon) * np.cos(lat)
    ds['Z'] = rad * np.sin(lat)

    ds['Vr'] = ds['V1'] * velocity_conversion
    ds['Density'] = ds['Density'] * density_conversion / \
        (1 - ds.attrs['xalpha']) * rad**2
    # Ram pressure (rho * v**2)
    ds['Pressure'] = ds['Density'] * ds['Vr']**2

    ds['X'].attrs.update({'units': 'AU'})
    ds['Y'].attrs.update({'units': 'AU'})
    ds['Z'].attrs.update({'units': 'AU'})
    ds = ds.rename_dims({'nevo': 'time'})
    ds = ds.assign_coords(time=t.data)

    # Remove unnecessary variables to make file sizes smaller
    ds = ds.drop_vars(['DT', 'NSTEP', 'BP', 'TIME',
                       'X1', 'X2', 'X3',
                       'B1', 'B2', 'B3',
                       'V1', 'V2', 'V3'])

    return ds


def process_directory(path):
    # Start off with evolution files
    fnames = sorted(glob.glob(path + '/evo.*.nc'))
    if len(fnames) == 0:
        raise ValueError("No evolution files found to process in the "
                         "current directory.")

    print("Beginning processing, may take several minutes")
    t0 = time.time()
    # Load and process the evo file
    datasets = [process_evo(xr.open_dataset(fname, engine='netcdf4'))
                for fname in fnames]

    # New path
    newpath = path + '/processed'
    if not os.path.exists(newpath):
        os.mkdir(newpath)

    for ds, fname in zip(datasets, fnames):
        # Save a single evolution file to NetCDF
        newfile = f"{newpath}/{os.path.basename(fname)}"
        ds.to_netcdf(newfile,
                     encoding={'time': {'units':
                                        'seconds since 1970-01-01'}})
        # Reload that file and save it to json
        # Loading with decode_times=False makes it so the datetimes
        # get encoded as ints within the JSON
        newds = xr.load_dataset(newfile, decode_times=False)

        # Drop the ".nc" extension and replace with json
        with open(newfile.replace('.nc', '.json'), 'w') as f:
            # Put it through dump/load/dump to get floating point truncation
            f.write(json.dumps(json.loads(json.dumps(newds.to_dict()),
                    parse_float=lambda x: round(float(x), 3))))

    # ---------
    # TIM file processing
    # ---------
    fnames = sorted(glob.glob(path + '/tim.*.nc'))
    if len(fnames) == 0:
        raise ValueError("No files found to process in the current directory.")

    print("Beginning processing, may take several minutes")
    t0 = time.time()
    # Load into a multi-file dataset to be able to concatenate along Time
    ds = xr.open_mfdataset(fnames, concat_dim='time', combine='by_coords',
                           preprocess=process_tim, engine='netcdf4')
    print(f"Dataset loaded: {time.time()-t0} s")

    ds.to_netcdf(f"{newpath}/test.nc", engine='scipy',
                 encoding={'time': {'units': 'seconds since 1970-01-01'}})
    # _, datasets = zip(*ds.groupby('time'))
    # paths = [f"{newpath}/tim.{i:04d}.nc" for i in range(len(fnames))]
    # xr.save_mfdataset(datasets, paths, engine='scipy')
    print(f"Dataset saved: {time.time()-t0} s")
    # for fname in fnames:
    #     with xr.open_dataset(fname) as ds:
    #         print(f"processing {fname}")
    #         process_ds(ds).to_netcdf(f'{newpath}/{os.path.basename(fname)}')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError("This program takes one argument, the path to the "
                         "directory to convert")
    process_directory(sys.argv[1])
