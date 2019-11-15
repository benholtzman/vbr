# fetches Q models from Iris.
#!pip install xarray

try:
    import urllib.request as urlrequest
except ImportError:
    import urllib as urlrequest
import xarray as xr # for loading netcdf
import os
import scipy.io as scp
import numpy as np

url_base='https://ds.iris.edu/files/products/emc/emc-files/'
iris_files={
 'Gung_Romanowicz_2002':
   {
	'server_name':'QRLW8_percent.nc',
	'dQ_field':'dq','z_field':'depth','lat_field':'latitude',
	'lon_field':'longitude','dims':'z,lat,lon'
   }
}

for ref in iris_files.keys():
    full_url=url_base+iris_files[ref]['server_name']
    if os.path.isfile(ref+'.nc') or os.path.isfile(ref+'.mat'):
        print(ref+' already downloaded.')
    else:
        print("attempting to fetch "+full_url)
        urlrequest.urlretrieve(full_url, ref+'.nc')
        print("file downloaded as ./"+ref+'.nc')


# slightly different fieldnames
for fi in iris_files.keys():
    if os.path.isfile(fi+'.mat') is False:
        ds=xr.open_dataset(fi+'.nc')
        
        if fi is 'Gung_Romanowicz_2002':
            print('ref Q is QL6c.1D')
            QL6c = np.tile(
                    np.array([[[70., 70, 70., 70., 80., 90., 100., 120., 130.,  
                                140., 150., 165., 165., 165., 165., 165., 165.,
                                165.]]]),
                    (91, 180, 1))
            Q_field = (
                ds[iris_files[fi]['dQ_field']].values.transpose(1, 2, 0)
                / 100 * QL6c + QL6c)
        else:
            Q_field = ds[iris_files[fi]['Q_field']].values.transpose(1, 2, 0)
        
        save_dict={'Latitude':ds[iris_files[fi]['lat_field']].values,
                   'Longitude':ds[iris_files[fi]['lon_field']].values,
                   'Depth':ds[iris_files[fi]['z_field']].values,
                   'Q':Q_field, 'Qinv':1/Q_field}
        print(fi+'.nc converted to '+fi+'.mat')
        scp.savemat(fi+'.mat',{'Q_Model':save_dict})
    else:
        print(fi+'.mat already exists')
