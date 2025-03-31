# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import metpy.calc as mpcalc
import cmocean
import cmcrameri
import cartopy.crs as ccrs
import os



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



# PARAMETERS
# Path
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'
data = '/lus/scratch/CT1/c1601279/rguillermin/GLORYS/raw_cli_mercator_Y2018_surf.nc'

# Selection
date = '2018-06-15'
SWIO = (25, 69, -36, 7)

# Plot
gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.5}
figsize = (20, 6)
velo_cmap = cmocean.cm.speed
vort_cmap = cmcrameri.cm.vik
ener_cmap = cmcrameri.cm.devon



d = xr.open_dataset(data)
u = d.uo
v = d.vo
time = d.time
d.close()
print('Data loaded.')

Lon, Lat = u.longitude, v.latitude
lon, lat = np.meshgrid(Lon, Lat)

dx = haversine(lon[:,:-1], lon[:,1:], lat[:,:-1], lat[:,1:])
dy = haversine(lon[:-1,:], lon[1:,:], lat[:-1,:], lat[1:,:])



if np.datetime64(date) not in time.values.astype('datetime64[D]'):
    print(f'{date} not found in simulation.')
    date = np.datetime_as_string(time.values.astype('datetime64[D]')[0])
    print(f'Setting {date} as the new date.')

years = np.unique(time.astype('datetime64[Y]')).astype('datetime64[Y]')
year = np.datetime64(date).astype('datetime64[Y]')


u = u.sel(time=str(year))
v = v.sel(time=str(year))



# Moyenne annuelle
u_yr = u.mean(dim='time')
v_yr = v.mean(dim='time')
print("Time mean calculated")

u = u.sel(time=date)
v = v.sel(time=date)
print("Data sliced")

# Vitesse turbulente
ut = (u - u_yr)
vt = (v - v_yr)
print("Turbulent velocity calculated")

u = u
v = v
    
u_yr = u_yr
v_yr = v_yr

velocity = np.sqrt(u ** 2 + v ** 2)
velocity_yr = np.sqrt(u_yr ** 2 + v_yr ** 2)
velocity_t = np.sqrt(ut ** 2 + vt ** 2)

# Calculate derivatives
dv_dlon = v.diff('longitude') / dx
du_dlat = u.diff('latitude') / dy
vorticity = dv_dlon - du_dlat

dv_yr_dlon = v_yr.diff('longitude') / dx
du_yr_dlat = u_yr.diff('latitude') / dy
vorticity_yr = dv_yr_dlon - du_yr_dlat

dvt_dlon = vt.diff('longitude') / dx
dut_dlat = ut.diff('latitude') / dy
vorticity_t = dvt_dlon - dut_dlat

KE = 1 / 2 * (u ** 2 + v ** 2)
MKE = 1 / 2 * (u_yr ** 2 + v_yr ** 2)
EKE = 1 / 2 * (ut ** 2 + vt ** 2)



fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(f"GLORYS {date}")

for ax in axs:
    ax.set_extent(SWIO)
    ax.coastlines(resolution='50m')
    ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black', zorder=3)
    ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5, zorder=3)
    ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5, zorder=3)

    land_color = ccrs.cartopy.feature.COLORS['land']

    minor_islands = ccrs.cartopy.feature.NaturalEarthFeature(
        category='physical',
        name='minor_islands',
        scale='10m',
        facecolor=land_color,
        edgecolor='black'
    )

    ax.add_feature(minor_islands, zorder=3)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {'color': 'k'}
    
    
# Velocity
ax = axs[0]
ax.set_title("Velocity")

a = -1.5
b = 0.2
norm = mpl.colors.LogNorm(10**a, 10**b)
data = velocity
cmap = velo_cmap
label = r'$v$ [m/s]'

pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')

# Vorticity
ax = axs[1]
ax.set_title("Vorticity")

norm = mpl.colors.Normalize(vmin=-5e-5, vmax=5e-5)
data = vorticity
cmap = vort_cmap
label = r'$\omega$ [h⁻¹]'
pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
ticks = cb.get_ticks()[1:-1]
cb.set_ticks(ticks)
cb.set_ticklabels([f'{tick:.0e}' for tick in ticks])

# Energy
ax = axs[2]
ax.set_title("Energy")

a = int(np.log10(np.nanmax(KE)))
b = a-2.5
norm = mpl.colors.LogNorm(vmin=10**b, vmax=10**a)
data = KE
cmap = ener_cmap
label = r'$KE$ [m².s⁻²]'

pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
    
fig.tight_layout()
filename = os.path.join(figures, f"dynamical_glorys_{date}.png")
fig.savefig(filename, dpi=300, transparent=True)
print(f'{filename} saved.')
plt.show()



fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(f"GLORYS")

for ax in axs:
    ax.set_extent(SWIO)
    ax.coastlines(resolution='50m')
    ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black', zorder=3)
    ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5, zorder=3)
    ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5, zorder=3)

    land_color = ccrs.cartopy.feature.COLORS['land']

    minor_islands = ccrs.cartopy.feature.NaturalEarthFeature(
        category='physical',
        name='minor_islands',
        scale='10m',
        facecolor=land_color,
        edgecolor='black'
    )

    ax.add_feature(minor_islands, zorder=3)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {'color': 'k'}
    
    
# Velocity
ax = axs[0]
ax.set_title("Velocity")

a = -1.5
b = 0.2
norm = mpl.colors.LogNorm(10**a, 10**b)
data = velocity_yr
cmap = velo_cmap
label = r'$v$ [m/s]'

pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')

# Vorticity
ax = axs[1]
ax.set_title("Vorticity")

norm = mpl.colors.Normalize(vmin=-5e-5, vmax=5e-5)
data = vorticity_yr
cmap = vort_cmap
label = r'$\omega$ [h⁻¹]'
pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
ticks = cb.get_ticks()[1:-1]
cb.set_ticks(ticks)
cb.set_ticklabels([f'{tick:.0e}' for tick in ticks])

# Energy
ax = axs[2]
ax.set_title("Energy")

a = int(np.log10(np.nanmax(MKE)))
b = a-2.5
norm = mpl.colors.LogNorm(vmin=10**b, vmax=10**a)
data = MKE
cmap = ener_cmap
label = r'$KE$ [m².s⁻²]'

pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
    
fig.tight_layout()
filename = os.path.join(figures, f"dynamical_glorys_{year}_mean.png")
fig.savefig(filename, dpi=300, transparent=True)
print(f'{filename} saved.')
plt.show()



fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(f"GLORYS Turbulent {date}")

for ax in axs:
    ax.set_extent(SWIO)
    ax.coastlines(resolution='50m')
    ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black', zorder=3)
    ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5, zorder=3)
    ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5, zorder=3)

    land_color = ccrs.cartopy.feature.COLORS['land']

    minor_islands = ccrs.cartopy.feature.NaturalEarthFeature(
        category='physical',
        name='minor_islands',
        scale='10m',
        facecolor=land_color,
        edgecolor='black'
    )

    ax.add_feature(minor_islands, zorder=3)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {'color': 'k'}
    
    
# Velocity
ax = axs[0]
ax.set_title("Velocity")

a = -1.5
b = 0.2
norm = mpl.colors.LogNorm(10**a, 10**b)
data = velocity_t
cmap = velo_cmap
label = r'$v$ [m/s]'

pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')

# Vorticity
ax = axs[1]
ax.set_title("Vorticity")

norm = mpl.colors.Normalize(vmin=-5e-5, vmax=5e-5)
data = vorticity_t
cmap = vort_cmap
label = r'$\omega$ [h⁻¹]'
pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
ticks = cb.get_ticks()[1:-1]
cb.set_ticks(ticks)
cb.set_ticklabels([f'{tick:.0e}' for tick in ticks])

# Energy
ax = axs[2]
ax.set_title("Energy")

a = int(np.log10(np.nanmax(EKE)))
b = a-2.5
norm = mpl.colors.LogNorm(vmin=10**b, vmax=10**a)
data = EKE
cmap = ener_cmap
label = r'$KE$ [m².s⁻²]'

pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
    
fig.tight_layout()
filename = os.path.join(figures, f"dynamical_glorys_{date}_turbulent.png")
fig.savefig(filename, dpi=300, transparent=True)
print(f'{filename} saved.')
plt.show()






