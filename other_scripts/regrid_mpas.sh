#!/bin/bash

export srcgrid=$1
export dstgrid=$2

export griddir=/glade/u/home/timand/remapping
export file=$3

export srcinitfile=${file}.nc
export dstinitfile=${file}.regrid.${dstgrid}.nc

echo "Extracting w, rho, theta"
ncks -v w,rho,theta ${srcinitfile} nCells.nc
echo "Renaming mislabeled dimensions"
ncrename -d nCells,ncol nCells.nc
echo "Extracing correctly dimensioned variables"
ncks -x -v w,rho,theta ${srcinitfile} ncol.nc
echo "Concatenating into output file."
ncks -A nCells.nc ncol.nc

ncremap -m ${griddir}/map_${srcgrid}_to_${dstgrid}_aave.nc -i\
 ncol.nc -o ${dstinitfile}