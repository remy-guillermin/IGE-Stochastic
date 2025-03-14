import xarray as xr
import numpy as np

ds = xr.open_dataset('/lus/work/CT1/c1601279/rguillermin/RUN_CROCO/run_swio2_deter2_2017_2023/swio_avg_pt1.nc')
u = ds.u[:,-1:,:,:]
v = ds.v[:,-1:,:,:]
w = ds.w[:,-1:,:,:]
salt = ds.salt[:,-1:,:,:]
temp = ds.temp[:,-1:,:,:]
zeta = ds.zeta
ds.close()

print('swio_avg_pt1.nc opened')

ds = xr.open_dataset('/lus/work/CT1/c1601279/rguillermin/RUN_CROCO/run_swio2_deter2_2017_2023/swio_avg_pt2.nc')
u2 = ds.u[:,-1:,:,:]
v2 = ds.v[:,-1:,:,:]
w2 = ds.w[:,-1:,:,:]
salt2 = ds.salt[:,-1:,:,:]
temp2 = ds.temp[:,-1:,:,:]
zeta2 = ds.zeta
ds.close()

print('swio_avg_pt2.nc opened')

ds = xr.open_dataset('/lus/work/CT1/c1601279/rguillermin/RUN_CROCO/run_swio2_deter2_2017_2023/swio_avg_pt3.nc')
u3 = ds.u[:,-1:,:,:]
v3 = ds.v[:,-1:,:,:]
w3 = ds.w[:,-1:,:,:]
salt3 = ds.salt[:,-1:,:,:]
temp3 = ds.temp[:,-1:,:,:]
zeta3 = ds.zeta
ds.close()

print('swio_avg_pt2.nc opened')

grid_path = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/CROCO_FILES/grid/croco_grid_swio2.nc'
g = xr.open_dataset(grid_path)
pm = g['pm'][:,:]
pn = g['pn'][:,:]
angle = g['angle'][:, :]
g.close()
print(f'{grid_path} opened')

u_tot = xr.concat([u, u2, u3], dim='time')
print('u concatenated')
v_tot = xr.concat([v, v2, v3], dim='time')
print('v concatenated')
w_tot = xr.concat([w, w2, w3], dim='time')
print('w concatenated')
salt_tot = xr.concat([salt, salt2, salt3], dim='time')
print('salinity concatenated')
temp_tot = xr.concat([temp, temp2, temp3], dim='time')
print('temperature concatenated')
zeta_tot = xr.concat([zeta, zeta2, zeta3], dim='time')
print('zeta concatenated')


angle_u = angle.rename({'eta_rho': 'eta_rho', 'xi_rho': 'xi_u'})
angle_u['xi_u'] = angle_u.xi_u + 0.5
angle_u['eta_rho'] = angle_u.eta_rho + 1

angle_v = angle.rename({'eta_rho': 'eta_v', 'xi_rho': 'xi_rho'})
angle_v['xi_rho'] = angle_v.xi_rho + 1
angle_v['eta_v'] = angle_v.eta_v + 0.5

u_cos = u_tot * xr.ufuncs.cos(angle_u)
u_cos['xi_u'] = u_cos.xi_u - 0.5
u_cos = u_cos.rename({'xi_u': 'xi_rho'})

v_sin = v_tot * xr.ufuncs.sin(angle_v)
v_sin['eta_v'] = v_sin.eta_v - 0.5
v_sin = v_sin.rename({'eta_v': 'eta_rho'})

u_geo = u_cos - v_sin
del u_cos, v_sin

u_sin = u_tot * xr.ufuncs.sin(angle_u)
u_sin['xi_u'] = u_sin.xi_u - 0.5
u_sin = u_sin.rename({'xi_u': 'xi_rho'})

v_cos = v_tot * xr.ufuncs.cos(angle_v)
v_cos['eta_v'] = v_cos.eta_v - 0.5
v_cos = v_cos.rename({'eta_v': 'eta_rho'})

v_geo = u_sin + v_cos

del u_sin, v_cos


u_geo = u_geo.reset_coords(drop=True)
v_geo = v_geo.reset_coords(drop=True)

u_geo.attrs = {
	'long_name': 'averaged geographic u-momentum component',
	'units': 'meter second-1',
	'field': 'u-geo-velocity, scalar, series',
	'standard_name': 'geographic_sea_water_x_velocity_at_rho_location'
}

v_geo.attrs = {
	'long_name': 'averaged geographic v-momentum component',
	'units': 'meter second-1',
	'field': 'v-geo-velocity, scalar, series',
	'standard_name': 'geographic_sea_water_y_velocity_at_rho_location'
}

velocity_h = np.sqrt(u_geo ** 2 + v_geo ** 2)
velocity_h.attrs = {
	'long_name': 'averaged geographic horizontal momentum component',
	'units': 'meter second-1',
	'field': 'geo-hor-velocity, scalar, series',
	'standard_name': 'geographic_sea_water_horizontal_velocity_at_rho_location'
	}


# Calculate derivatives
dv_dlon = v_geo.differentiate('xi_rho') * pm[:-1,:-1]
du_dlat = u_geo.differentiate('eta_rho') * pn[:-1,:-1]

vorticity = dv_dlon - du_dlat

vorticity.attrs = {
	'long_name': 'averaged geographic vorticity',
	'units': 'second-1',
	'field': 'vorticity, scalar, series',
	'standard_name': 'geographic_sea_water_vorticity_at_rho_location'
}


ds = xr.Dataset({
	'u': u_tot.reset_coords(drop=True),
	'u_geo': u_geo,
	'v': v_tot.reset_coords(drop=True),
	'v_geo': v_geo,
	'w': w_tot.reset_coords(drop=True),
	'V_hor': velocity_h,
	'vort': vorticity,
	'salt': salt_tot.reset_coords(drop=True),
	'temp': temp_tot.reset_coords(drop=True),
	'zeta': zeta_tot.reset_coords(drop=True)
})

ds.to_netcdf('/lus/work/CT1/c1601279/rguillermin/RUN_CROCO/run_swio2_deter2_2017_2023/swio_avg_suf.nc')
	
	
	