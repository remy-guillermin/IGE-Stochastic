import xarray as xr
import numpy as np
import glob

def haversine(lon1, lon2, lat1, lat2):
    lon1 = np.deg2rad(lon1)
    lon2 = np.deg2rad(lon2)
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.atan2(np.sqrt(a), np.sqrt(1-a))
    
    r = 6371000
    return c * r

files = glob.glob('/lus/work/CT1/c1601279/rguillermin/GLORYS/raw_cli_mercator_Y*.nc')
files = np.sort(files)

ds = xr.open_dataset(files[0])
u = ds.uo[:,0:1,:,:]
v = ds.vo[:,0:1,:,:]
ds.close()
print(f'Opened {files[0]}')

for file in files[1:]:
    ds = xr.open_dataset(file)
    u = xr.concat([u, ds.uo[:,0:1,:,:]], dim='time')
    print('u concatenated')
    v = xr.concat([v, ds.vo[:,0:1,:,:]], dim='time')
    print('v concatenated')
    ds.close()
    print(f'Opened {file}')
       
Lon = u.longitude
Lat = u.latitude

lon, lat = np.meshgrid(Lon, Lat)

dx = haversine(lon[:,:-1], lon[:,1:], lat[:,:-1], lat[:,1:])
dy = haversine(lon[:-1,:], lon[1:,:], lat[:-1,:], lat[1:,:])
print(dx.shape, dy.shape)

velocity_h = np.sqrt(u**2 + v**2)

velocity_h.attrs = {
    'cell_methods': 'area: mean',
	'long_name': 'Surface horizontal velocity',
	'units': 'm s-1',
 	'unit_long': 'Meters per second',
	'standard_name': 'horizontal_sea_water_velocity',
}


# Calculate derivatives
dv_dlon = v[:,:,:-1,:-1].differentiate('longitude') / dx[:-1,:]
du_dlat = u[:,:,:-1,:-1].differentiate('latitude') / dy[:,:-1]
vorticity = dv_dlon - du_dlat

vorticity.attrs = {
    'cell_methods': 'area: mean',
	'long_name': 'Surface horizontal vorticity',
	'units': 's-1',
	'unit_long': 'per second',
	'standard_name': 'sea_water_vorticity'
}

ds = xr.Dataset({
	'u': u,
	'v': v,
	'V_hor': velocity_h,
	'vort': vorticity,
})

ds.to_netcdf('/lus/work/CT1/c1601279/rguillermin/GLORYS/glorys_avg_suf.nc')

