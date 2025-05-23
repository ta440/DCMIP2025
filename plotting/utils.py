import matplotlib.colors as colors
import numpy as np

# A function to perform interpolation to a constant altitude

def z_interp(h, field_vals_all, lon, lat, z_val):
    field_vals = np.zeros((len(lat), len(lon)))
    for i in np.arange(len(lat)):
        for j in np.arange(len(lon)):
            if h[-1, i, j] > z_val:
                # This value is inside the topography
                # field_vals[i, j] = np.nan
                field_vals[i, j] = field_vals_all[-1,i,j]
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

# interpolates u and v at the same time
def z_interp_uv(h, u_field, v_field, lon, lat, z_val):
    field_vals = np.zeros((2, len(lat), len(lon)))
    for i in np.arange(len(lat)):
        for j in np.arange(len(lon)):
            if h[-1, i, j] > z_val:
                # This value is inside the topography
                # field_vals[0, i, j] = np.nan
                # field_vals[1, i, j] = np.nan
                field_vals[0, i, j] = u_field[-1, i, j]
                field_vals[1, i, j] = v_field[-1, i, j]
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


class MidpointNormalize(colors.Normalize):  # a class for the colorbar normalization
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        v_ext = np.max([np.abs(self.vmin), np.abs(self.vmax)])
        x, y = [-v_ext, self.midpoint, v_ext], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))