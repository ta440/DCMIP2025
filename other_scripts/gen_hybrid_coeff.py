# -*- coding: utf-8 -*-
"""
Created on Thu 7 Nov, 2024

Generate new hybrid coordinates to feed into the CAM model.
Here, we generate a vertical grid for a constant spacing
and fixed height assuming an isothermal atmosphere.

@author: ta440 and nandroski
"""

import numpy as np
import scipy
import csv
import matplotlib.pyplot as plt
from os.path import abspath, dirname

# Model parameters
z_top = 30000 # Model top (m)
dz_low = 300 # Constant vertical spacing (m)
dz_max = 1000
h0 = 1500 # Mountain height
p0 = 1000 # Reference pressure (hPa)
T0 = 288 # Isothermal temperature (K)
g = 9.80616 # Gravitational acceleration (ms^{-2})
Rd = 287.04  # Universal gas constant

# save file parameters
modname = ''
desc = ''
attribution = 'DCMIP 2025'
creator = 'Nicholas Androski'
output_dir = ''

print(f'Lowest dz is {dz_low}')

# Guess a geometric stretching factor:
# 1.02 gets dz = 100 to 1000 at z = 10km
# 1.009 gets dz = 100 to 500 at z = 10km
fact = 1.02

#-------------------------------------------------
# Step 1: Compute z values
#-------------------------------------------------

# Region 1:
# We want increasing dz from 100 to 800 for the first 21 km,
# then keep at 800 for the remainder of the domain.
z_vals = [0]
dz_vals = []

while z_vals[-1] < h0:
    z_vals.append(z_vals[-1]+dz_low)
    dz_vals.append(dz_low)

print(f'There are {len(dz_vals)} levels over the mountain')

# Region 2:
# For dz=100, it will stretch to hit dz_max
# at 10 km.
while dz_vals[-1] < dz_max:
    dz_new = dz_vals[-1]**fact
    if dz_new > dz_max:
        dz_new = dz_max
    
    dz_vals.append(dz_new)
    z_vals.append(z_vals[-1]+dz_new)

print(f'We have stopped stretching at z = {z_vals[-1]} after {len(dz_vals)} levels')

# Region 3:
# Increase up to z = 30000
while z_vals[-1] < z_top:
    z_vals.append(z_vals[-1]+dz_max)
    dz_vals.append(dz_max)


#print(f'dz vals are {dz_vals}')
#print(f'z vals are {z_vals}')


print(f'There are {len(z_vals)} interface levels')
print(f'There are {len(dz_vals)} midpoint levels')
print(f'z top is {z_vals[-1]}')
#print(z_vals)

z_vals = np.asarray(z_vals)
print(z_vals)

#-------------------------------------------------
# Step 2: Compute p values at the interfaces.
#-------------------------------------------------

# Isothermal p-z relation
def p_isothermal(z):
    return p0*np.exp((-g*z)/(Rd*T0))

pi_vals = p_isothermal(z_vals)

print('p values are')
print(pi_vals)


#-------------------------------------------------
# Step 3: Compute A and B coefficients 
#-------------------------------------------------

# Transition exponent between p0 and ps
# c=1 resembles the sigma-coordinate system
# c=2 allows a more gradual transition from
# contour-following to flat pressure levels
c = 1

# eta coordinate
eta = pi_vals/p0

# Compute interface values. These are equations
# 124 and 125 from DCMIP2008 and Laprise, Girard 1990.
hybi = ((eta - eta[-1])/(1-eta[-1]))**c
hyai = eta - hybi

# Compute midpoint values
hyam = np.zeros(len(hyai) - 1)
hybm = np.zeros(len(hybi) - 1)

for i in np.arange(len(hyam)):
    hyam[i] = 0.5*(hyai[i+1]+hyai[i])
    hybm[i] = 0.5*(hybi[i+1]+hybi[i])

# Compute the pressures at the midpoints
pm_vals = np.zeros_like(hyam)

ps = 1000.
for i in np.arange(len(pm_vals)):
    pm_vals[i] = hyam[i]*p0 + hybm[i]*ps

# Check that we have agreement of pressures at the
# interfaces
diff = 0
for i in np.arange(len(pi_vals)):
    pi_check = hyai[i]*p0 + hybi[i]*ps
    diff += np.abs(pi_check-pi_vals[i])

print(f'difference between computed pi values is {diff}')

#-------------------------------------------------
# Step 4: Save values in a netCDF file
#-------------------------------------------------

save = True

filename = f'{abspath(dirname(__file__))}/coord_text_files/mount_horiz_stretched_coords_T0_{T0}_dz_low_{dz_low}_dz_max_{dz_max}_ztop_{z_top}_c_{c}.txt'
delimiter = ', '

# Values need to be saved in order of increasing pressure
pm_vals = np.flip(pm_vals)
hyam = np.flip(hyam)
hybm = np.flip(hybm)
pi_vals = np.flip(pi_vals)
hyai = np.flip(hyai)
hybi = np.flip(hybi)

ilev = pi_vals/p0 * 1000
lev = pm_vals/p0 * 1000


def create_vertGrid_ncfile(ilev, lev, hyai, hybi, hyam, hybm, P0, PS,
                           output_dir='', modname='',desc='',
                           attribution='DCMIP 2025', creator = 'Nicholas Androski'):
    from netCDF4 import Dataset
    nlev = len(lev)
    
    ncfile = Dataset(f"{output_dir}template_L{int(nlev)}_{modname}.nc",
                     mode='w',data_model="NETCDF3_CLASSIC")
    
    # Define dimensions
    lev_dim = ncfile.createDimension("lev", nlev)      # Fixed size
    ilev_dim = ncfile.createDimension("ilev", nlev+1)     # Fixed size
    
    # Create coordinate variables
    levs = ncfile.createVariable("lev", np.float32, ("lev",))
    ilevs = ncfile.createVariable("ilev", np.float32, ("ilev",))

    # reference variables
    P0_ = ncfile.createVariable("P0", np.float64)  # Scalar
    PS_ = ncfile.createVariable("PS", np.float64, ("lev",))  # Surface pressure field
    
    # Create a data variable
    hyai_ = ncfile.createVariable("hyai", np.float32, ("ilev"), 
                                 zlib=True, fill_value=-9999.0)
    hybi_ = ncfile.createVariable("hybi", np.float32, ("ilev"), 
                                 zlib=True, fill_value=-9999.0)
    hyam_ = ncfile.createVariable("hyam", np.float32, ("lev"), 
                                 zlib=True, fill_value=-9999.0)
    hybm_ = ncfile.createVariable("hybm", np.float32, ("lev"), 
                                 zlib=True, fill_value=-9999.0)

    levs.A_var = "hyam"
    levs.B_var = "hybm"
    levs.P0_var = "P0"
    levs.PS_var = "PS"
    levs.bounds = "ilev"

    ilevs.A_var = "hyai"
    ilevs.B_var = "hybi"
    ilevs.P0_var = "P0"
    ilevs.PS_var = "PS"
    
    # Add global attributes following NCAR conventions
    ncfile.title = f"L{nlev} Hybrid eta coefficients with spacing {modname}"
    ncfile.institution = attribution
    ncfile.source = f"Derived using {desc}"
    ncfile.history = "Created " + Dataset.__dict__.get('__name__', creator) # add your name here
    
    # Add variable attributes
    levs.units = "hybrid_sigma_pressure (hPa)"
    levs.long_name = "hybrid level at layer midpoints (1000*(A+B))" 
    
    levs.units = "hybrid_sigma_pressure (hPa)"
    levs.long_name = "hybrid level at layer interfaces (1000*(A+B))" 
    
   
    hyai_.long_name = "hybrid A coefficient at layer interfaces"
    hybi_.long_name = "hybrid B coefficient at layer interfaces"
    hyam_.long_name = "hybrid A coefficient at layer midpoints"
    hybm_.long_name = "hybrid B coefficient at layer midpoints"
    
    # Populate coordinate variables
    if lev[0]>lev[-1]:
        levs[:] = lev[::-1] #flip all the arrays so that levels go from low to high
        ilevs[:] = ilev[::-1]
    else:
        levs[:] = lev #flip all the arrays so that levels go from low to high
        ilevs[:] = ilev
        
        

    # Populate ref variables
    P0_.assignValue(P0) 
    
    # Populate coefficient variables
    if lev[0]>lev[-1]:
        hyai_[::] = hyai[::-1]
        hybi_[::] = hybi[::-1]
        hyam_[::] = hyam[::-1]
        hybm_[::] = hybm[::-1]
    else:
        hyai_[:] = hyai
        hybi_[:] = hybi
        hyam_[:] = hyam
        hybm_[:] = hybm
        
    
    # Close the NetCDF file
    ncfile.close()

    return

if save == True:
    create_vertGrid_ncfile(ilev, lev, hyai, hybi, hyam, hybm, p0, ps,
                  output_dir=output_dir, modname=modname, desc=desc,
                 attribution=attribution, creator=creator)
    
    print('saved output to', f"{output_dir}template_L{len(lev)}_{modname}.nc")





