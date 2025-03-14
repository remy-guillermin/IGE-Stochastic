import xarray as xr
import numpy as np

def haversine(lon1, lon2, lat1, lat2):
    lon1 = np.deg2rad(lon1)
    lon2 = np.deg2rad(lon2)
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.atan2(np.sqrt(a), np.sqrt(1-a))
    
    r = 6371
    return c * r

ds = xr.open_dataset('/home/guilremy/IGE-Stochastic/copernicus-data/CMEMS_Indian_currents_surface.nc')

u, v = ds.uo, ds.vo

ds.close()

Lon = u.longitude
Lat = u.latitude

lon, lat = np.meshgrid(Lon, Lat)

dx = haversine(lon[:,:-1], lon[:,1:], lat[:,:-1], lat[:,1:])
dy = haversine(lon[:-1,:], lon[1:,:], lat[:-1,:], lat[1:,:])
print(dx.shape, dy.shape)