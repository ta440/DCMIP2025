# DCMIP2025
A collection of scripts for the DCMIP2025 test cases. There are three test cases:
1. Gravity wave breaking over two mountains
2. Mountain-driven mesoscale phenomena: Gap flow and vortex shedding
3. Squall line with Kessler physics and radar reflectivity

Here is an overview of each directory:

## CAM_src_code
New src code for the tests, inluding the initial condition scripts, kessler physics routines.

## RF_SourceMods
SourceMods for the Rayleigh friction damping layers, which were implemented using the FHS94 compset.

## other_scripts
Random stuff, such as regridding files.

## user_nl_cam
Templates for user_nl_cam for each test and dynamical core combination.

## xmlchanges
The xmlchanges used to set up the tests.

