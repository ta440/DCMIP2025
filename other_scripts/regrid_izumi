#!/bin/bash

export srcgrid=$1
export dstgrid=$2

export griddir=/project/amp02/dcmip2025/weight_files
export file=$3

export srcinitfile=${file}.nc
export dstinitfile=${file}.regrid.${dstgrid}.nc

ncremap -m ${griddir}/map_${srcgrid}_to_${dstgrid}_aave.nc -i\
 ${srcinitfile} -o ${dstinitfile}
