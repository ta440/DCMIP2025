#!/bin/bash

export srcgrid=$1
export dstgrid=$2

export griddir=/glade/u/home/timand/remapping
export file=$3

# Set up a temporary space for it so we avoid file system contention
export tmpspace=$(mktemp -d -p /dev/shm/ )
mkdir -p $tmpspace
echo "Using temp space: ${tmpspace}"

export srcinitfile=${file}.nc
export dstinitfile=${file}.regrid.${dstgrid}.nc

echo "Extracting w, rho, theta"
time ncks -v w,rho,theta ${srcinitfile} ${tmpspace}/nCells.nc
echo "Renaming mislabeled dimensions"
time ncrename -d nCells,ncol ${tmpspace}/nCells.nc
echo "Extracting correctly dimensioned variables"
time ncks -x -v w,rho,theta ${srcinitfile} ${tmpspace}/ncol.nc
echo "Concatenating into output file."
time ncks -A ${tmpspace}/nCells.nc ${tmpspace}/ncol.nc
time ncremap -m ${griddir}/map_${srcgrid}_to_${dstgrid}_aave.nc -i ${tmpspace}/ncol.nc -o ${dstinitfile}

# Clean up the temporary space:
rm -rf ${tmpspace}