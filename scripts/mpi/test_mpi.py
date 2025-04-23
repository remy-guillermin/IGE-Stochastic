from mpi4py import MPI
import xarray as xr
import numpy as np
import glob
import os

def calc_depth(
    s: xr.DataArray, 
    Cs: xr.DataArray, 
    hc: float, 
    h: xr.DataArray
) -> xr.DataArray:
    """
    Compute depth using the S-coordinate transformation.

    Parameters
    ----------
    s : xr.DataArray
        S-coordinate at RHO-points, typically ranging from -1 (surface) to 0 (bottom).
    Cs : xr.DataArray
        S-coordinate stretching curves at W-points, defining the vertical stretching.
    hc : float
        Critical depth parameter (in meters), influencing vertical terrain-following transformation.
    h : xr.DataArray
        Bathymetric depth at RHO-points (in meters), representing the seafloor depth.

    Returns
    -------
    xr.DataArray
        Computed depth at RHO-points.
    """
    N = len(s)
    M, L = h.shape
    z0 = np.zeros((N, M, L))
    depth = np.zeros((N, M, L))
    for k in range(N):
        z0[k, :, :] = (hc * s[k] + h * Cs[k]) / (hc + h)
        depth[k, :, :] = z0[k, :, :] * h
    return depth

CRT1 = 0.03 # kg/m3
    
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

scratch = '/lus/scratch/CT1/c1601279/rguillermin'
scratch_lisa = '/lus/scratch/CT1/c1601279/lweiss'
work = '/lus/work/CT1/c1601279/rguillermin'
store = '/lus/store/CT1/c1601279/rguillermin'
stoens = 'run_swio2_stoens30_2018_CD'
grid = '/lus/work/CT1/c1601279/rguillermin/grid/croco_grid_swio2.nc'

g = xr.open_dataset(grid)[['lon_rho', 'lat_rho', 'mask_rho', 'h']]
lon = g.lon_rho
lat = g.lat_rho
eta_rho = g.eta_rho
xi_rho = g.xi_rho
h = g.h
mask_rho = g.mask_rho
g.close()

ensemble_path = sorted(glob.glob(os.path.join(scratch_lisa, 'CROCO', 'SWIOSE_dev', stoens, '*avg.nc')))
if rank == 0:
    print(f'{len(ensemble_path)} files found for {stoens}')

ds = xr.open_dataset(ensemble_path[0])[['s_rho', 'Cs_rho', 'hc']]
depth_sigma = calc_depth(ds.s_rho, ds.Cs_rho, ds.hc, h)
s_rho = ds.s_rho
ds.close()

depth_level = xr.Dataset(
    data_vars=dict(
        depth_level=(['s_rho', 'eta_rho', 'xi_rho'], depth_sigma),
    ),
    coords=dict(
        eta_rho=eta_rho,
        xi_rho=xi_rho,
        s_rho=s_rho,
    )
)

h = h.where((h != 50), 0)

files_per_proc = np.array_split(ensemble_path, size)

my_files = files_per_proc[rank]

results = []

for path in my_files:
    print(f'Rank {rank} processing {path}')
    
    ds = xr.open_dataset(path)['rho']
    

    MLDC = ds - ds.isel(s_rho=-2) < CRT1 # Condition on density for MLD
    ds.close()
    mask = - MLDC.sum(dim='s_rho') 
    mask['eta_rho'] = mask.eta_rho - 1
    mask['xi_rho'] = mask.xi_rho - 1
    MLD = depth_level.isel(s_rho=mask)
    MLD = MLD.where((g.mask_rho), np.nan)
    
    results.append(MLD)
    
    
if results:
    local_result = xr.concat(results, dim='member')
else:
    local_result = None
    
gathered = comm.gather(local_result, root=0)

if rank==0:
    final_results = xr.concat([ds for ds in gathered if ds is not None], dim='member')
    
    print(final_results)