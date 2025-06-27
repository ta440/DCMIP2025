import matplotlib.colors as colors
import numpy as np
import xarray as xr

# A function to perform interpolation to a constant altitude

def z_interp(h, field_vals_all, lon, lat, z_val):
    field_vals = np.zeros((len(lat), len(lon)))
    for i in np.arange(len(lat)):
        for j in np.arange(len(lon)):
            if h[-1, i, j] > z_val:
                # This value is inside the topography
                field_vals[i, j] = np.nan
                # field_vals[i, j] = field_vals_all[-1,i,j]
            else:
                # Find indices either side of this value
                low_idx = np.where(h[:, i, j] < z_val)[0][0]
                high_idx = np.where(h[:, i, j] > z_val)[0][-1]

                # Compute weightings
                weight_low = (z_val - h[low_idx, i, j])/(h[high_idx, i, j] - h[low_idx, i, j])
                weight_high = 1. - weight_low
    
                # Compute and store value
                field_vals[i, j] = weight_low*field_vals_all[low_idx, i, j] + weight_high * field_vals_all[high_idx, i, j]
    return field_vals

# interpolates rad_ref for an array of z_vals, can be used to compute the base reflectivity or the reflectivity for a radar scan
def z_interp_br(h, rad_ref_field, z_val):
    # h and rad_ref are given with flattened latlon arrays
    field_vals = np.zeros(rad_ref_field.shape[1])
    for i in np.arange(len(field_vals)):
            if (h[-1,i] > z_val[i]):
                # this value is below the lowest model level/topography and would require extrapolation
                field_vals[i] = np.nan
                
            else:
                # Find indices either side of this value
                low_idx = np.where(h[:, i] < z_val[i])[0][0]
                high_idx = np.where(h[:, i] > z_val[i])[0][-1]
                
                # Compute weightings
                weight_low = (z_val[i] - h[low_idx, i])/(h[high_idx, i] - h[low_idx, i])
                weight_high = 1. - weight_low
    
                # Compute and store value
                field_vals[i] = weight_low*rad_ref_field[low_idx, i] + weight_high * rad_ref_field[high_idx, i]
    return field_vals

# interpolates u and v at the same time
def z_interp_uv(h, u_field, v_field, lon, lat, z_val):
    field_vals = np.zeros((2, len(lat), len(lon)))
    for i in np.arange(len(lat)):
        for j in np.arange(len(lon)):
            if h[-1, i, j] > z_val:
                # This value is inside the topography
                field_vals[0, i, j] = np.nan
                field_vals[1, i, j] = np.nan
                # field_vals[0, i, j] = u_field[-1, i, j]
                # field_vals[1, i, j] = v_field[-1, i, j]
            else:
                # Find indices either side of this value
                low_idx = np.where(h[:, i, j] < z_val)[0][0]
                high_idx = np.where(h[:, i, j] > z_val)[0][-1]

                # Compute weightings
                weight_low = (z_val - h[low_idx, i, j])/(h[high_idx, i, j] - h[low_idx, i, j])
                weight_high = 1. - weight_low

                # Compute and store value
                field_vals[0, i, j] = weight_low*u_field[low_idx, i, j] + weight_high * u_field[high_idx, i, j]
                field_vals[1, i, j] = weight_low*v_field[low_idx, i, j] + weight_high * v_field[high_idx, i, j]
    return field_vals

def z_interp_w_hydrostatic(nc, time, lon, lat, z_val):
    h = nc['Z3'][time, :, lat, lon]
    try:
        q = nc['Q'][time, :, lat, lon]
    except:
        q = 0
        qval = 0
    T_field = nc['T'][time, :, lat, lon]
    hyam = nc['hyam']
    hybm = nc['hybm']
    PS = nc['PS'][time, lat, lon]
    omega_field = nc['OMEGA'][time, :, lat, lon]
    field_vals = np.zeros((len(lat), len(lon)))
    for i in np.arange(len(lat)):
        for j in np.arange(len(lon)):
            if h[-1, i, j] > z_val:
                # This value is inside the topography
                field_vals[i, j] = np.nan
            else:
                # Find indices either side of this value
                low_idx = np.where(h[:, i, j] < z_val)[0][0]
                high_idx = np.where(h[:, i, j] > z_val)[0][-1]

                # Compute weightings
                weight_low = (z_val - h[low_idx, i, j])/(h[high_idx, i, j] - h[low_idx, i, j])
                weight_high = 1. - weight_low

                # Compute and store value
                if (q != 0):
                    qval = weight_low * q[low_idx, i, j] + weight_high * q[high_idx, i, j]
                T = weight_low * T_field[low_idx, i, j] + weight_high * T_field[high_idx, i, j]
                p1 = hyam[low_idx]*1e5 + hybm[low_idx]*PS[i, j]
                p2 = hyam[high_idx]*1e5 + hybm[high_idx]*PS[i, j]
                P = weight_low * p1 + weight_high * p2
                Tv = T * (1.0 + 0.608 * qval)
                rho = P / (287.0 * Tv)
                omega = weight_low * omega_field[low_idx, i, j] + weight_high * omega_field[high_idx, i, j]
                field_vals[i, j] = -omega/(rho * 9.81)
    return field_vals

def z_interp_w_nonhydro(nc, time, lon, lat, z_val):
    h = nc['Z3'][time, :, lat, lon]
    try:
        q = nc['Q'][time, :, lat, lon]
        qval = 1
    except:
        q = 0
        qval = 0
    T_field = nc['T'][time, :, lat, lon]
    hyam = nc['hyam']
    hybm = nc['hybm']
    P_field = nc['PMID'][time, :, lat, lon]
    omega_field = nc['OMEGA'][time, :, lat, lon]
    field_vals = np.zeros((len(lat), len(lon)))
    for i in np.arange(len(lat)):
        for j in np.arange(len(lon)):
            if h[-1, i, j] > z_val:
                # This value is inside the topography
                field_vals[i, j] = np.nan
            else:
                # Find indices either side of this value
                low_idx = np.where(h[:, i, j] < z_val)[0][0]
                high_idx = np.where(h[:, i, j] > z_val)[0][-1]

                # Compute weightings
                weight_low = (z_val - h[low_idx, i, j])/(h[high_idx, i, j] - h[low_idx, i, j])
                weight_high = 1. - weight_low

                # Compute and store value
                if (qval != 0):
                    qval = weight_low * q[low_idx, i, j] + weight_high * q[high_idx, i, j]
                T = weight_low * T_field[low_idx, i, j] + weight_high * T_field[high_idx, i, j]
                P = weight_low * P_field[low_idx, i, j] + weight_high * P_field[high_idx, i, j]
                Tv = T * (1.0 + 0.608 * qval)
                rho = P / (287.0 * Tv)
                omega = weight_low * omega_field[low_idx, i, j] + weight_high * omega_field[high_idx, i, j]
                field_vals[i, j] = -omega/(rho * 9.81)
    return field_vals

def z_interp_xr(ds, t_idxs, z_levs, lon_inds):
    out = np.zeros((len(t_idxs), len(z_levs), len(lon_inds)))
    for i in range(len(t_idxs)):
        for j in range(len(lon_inds)):
            slice = ds.isel(time=t_idxs[i], lon=lon_inds[j])
            out[i, :, j] = slice.set_index(lev='Z3').interp(lev=z_levs).to_dataarray().data
    return out

def z_interp_divvor(divvor, t_idxs, z_levs, lon_inds):
    out = np.zeros((len(t_idxs), len(z_levs), len(lon_inds)))
    for i in range(len(t_idxs)):
        for j in range(len(lon_inds)):
            slice = divvor.isel(time=i, lon=j)
            out[i, :, j] = slice.set_index(lev='Z3').interp(lev=z_levs).to_dataarray().data
    return out


class MidpointNormalize(colors.Normalize):  # a class for the colorbar normalization
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        v_ext = np.max([np.abs(self.vmin), np.abs(self.vmax)])
        x, y = [-v_ext, self.midpoint, v_ext], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def base_reflectivity(nc, t_idxs, lat_idxs, lon_idxs,
                      h0, lat0, lon0, epsilon0, max_range,
                      R_earth, small_Earth_factor, fullEarthRadar_smallEarthDistance=False):
    '''
    Calculates the base reflectivity, the radar reflectivity of the 
    lowest elevation (epsilon0) radar scan from a weather radar located
    at (lat0, lon0) at a height h0 above Earth's surface with the maximum
    along beam radar range given by max_range.
    Uses the 43ERM to approximate the radar beam path, which assumes an
    effective curvature of 4/3 R_earth of the beam from Earth's surface.
    ------------------------------------
    Inputs:
    ------------------------------------
    nc                 : The netCDF dataset containing lat, lon, RAD_REF, and Z3. Assumes lat-lon grid.
    
    t_idxs             : time indices to calculate base reflectivity
    
    lat_idxs           : latitude indices to calculate base reflectivity
    
    lon_idxs           : longitude indices to calculate base reflectivity
    
    h0                 : The height of the simulated radar above the Earth's surface in meters
    
    lat0               : The latitude location of the simulated radar, in degrees
    
    lon0               : The longitude location of the simulated radar, in degrees
    
    epsilon0           : The elevation (angle above Earth's local surface) of the radar scan in degrees
    
    max_range          : The maximum along beam distance that the radar beam travels in meters.
    
    R_earth            : The radius of the actual Earth in meters (6.371220E6 recommended)
    
    small_Earth_factor : The reduced Earth scale factor such that the radius is given by R_earth/small_Earth_factor
    
    fullEarthRadar_smallEarthDistance: If True, uses R_earth/small_Earth_factor to calculate 
                                       great circle distances between (lat, lon) points but R_earth
                                       to calculate the curvature of the radar beam.
    ------------------------------------
    Output:
    ------------------------------------
    br : The beam reflectivity of shape (len(t_idxs),len(lat_idxs),len(lon_idxs))

    '''

    if fullEarthRadar_smallEarthDistance:
        R_earth_for_distances = R_earth/small_Earth_factor
        R_earth_for_beam      = R_earth
    else:
        R_earth_for_distances = R_earth/small_Earth_factor
        R_earth_for_beam      = R_earth/small_Earth_factor

    
    def great_circle_distance(lat1, lon1, lat2, lon2):
        lat1, lon1, lat2, lon2 = np.deg2rad(lat1), np.deg2rad(lon1), np.deg2rad(lat2), np.deg2rad(lon2)
        dlon = lon2 - lon1
        d = R_earth_for_distances * np.arccos(np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(dlon))
        return d
    
    def compute_s_from_latlon(lat, lon, lat0, lon0):
        s = great_circle_distance(lat0, lon0, lat, lon)
        return s
    
    def radar_path(h0, epsilon0, r):
        Reff = 4/3 * R_earth_for_beam
        epsilon0_rad = np.radians(epsilon0)
        h = np.sqrt((Reff + h0)**2 + r**2 + 2*(Reff+h0)*r*np.sin(epsilon0_rad)) - Reff
        s = Reff * np.arcsin(r*np.cos(epsilon0_rad) / (Reff + h))
        
        return h, s

    def find_h_from_s(s):
        from scipy.optimize import root_scalar
        
        def func_to_solve(r):
            h_test, s_test = radar_path(h0, epsilon0, r)
            return s - s_test
            
        try:
            sol = root_scalar(func_to_solve, bracket=[0, max_range], method='brenth',xtol=eps,maxiter=iter_num)
            if sol.converged:
                h_sol, s_sol = radar_path(h0, epsilon0, sol.root)
                if verbose:
                  print(f"The h value corresponding to s = {s} is approximately h = {h_sol}")
            else:
                print("Root finding did not converge.")
        except ValueError as e:
            print(f"Error encountered: {e}")

        return h_sol

    
    # setup variables
    rad_ref  = nc['RAD_REF']
    Z3       = nc['Z3']
    lats     = nc['lat'][lat_idxs]
    lons     = nc['lon'][lon_idxs]

    # settings for h_from_s iteration
    h_from_s_vec = np.vectorize(find_h_from_s)
    eps = 1E-12
    iter_num = 100
    verbose = False

    # initialize base reflectivity
    br = np.ones((len(t_idxs), len(lats), len(lons))) * np.nan

    # compute the maximum height and distances sampled by the radar pulse
    h_max, s_max = radar_path(h0, epsilon0, max_range)
    print(f'h_max: {h_max/1.0e3:.2f} km')
    print(f's_max: {s_max/1.0e3:.2f} km')

    
    k = 0
    for t_idx in t_idxs:
        print(f'Time Index {t_idx}')
        
        # compute the along great circle distances between (lat, lon) and (lat0, lon0)
        s = np.ones((len(lats), len(lons)))
        for i in range(0,len(lats)):
            for j in range(0, len(lons)):
                s[i,j] = compute_s_from_latlon(lats[i], lons[j], lat0, lon0)
        
        # compute in range and out of range masks for s
        # if s is larger than s_max then it is out of range of the radar
        s_mask = np.where(s > s_max, True, False)
        # both s_mask and h_mask should be of shape (lat, lon)
        out_of_range_mask = s_mask
        in_range_mask = np.logical_not(out_of_range_mask)

        # compute the height of the radar beam path for lat, lon points that are in range
        h = np.ones((len(lats), len(lons))) * np.nan
        h[in_range_mask] = h_from_s_vec(s[in_range_mask])

        # compute a combined out of range mask accounting for s and h
        # any columns in which h is above z_top are out of range
        h_mask = np.where((np.max(Z3[t_idx, :, lat_idxs, lon_idxs],axis=0) < h) | (h == np.nan), True, False)
        
        out_of_range_mask = h_mask | s_mask
        in_range_mask = np.logical_not(out_of_range_mask)

        # select in range Z3 and rad_ref
        Z3_mask = Z3[t_idx, :, lat_idxs, lon_idxs]
        Z3_mask = Z3_mask[:, in_range_mask]
        rad_ref_mask = rad_ref[t_idx, :, lat_idxs, lon_idxs]
        rad_ref_mask = rad_ref_mask[:, in_range_mask]

        # interpolate rad_ref to the given h for each in range (lat, lon)
        br[k, in_range_mask] = z_interp_br(Z3_mask, rad_ref_mask, h[in_range_mask])
        k+=1
    return br