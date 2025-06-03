for i in 88 120 207
do
echo "Extracting w, rho, theta for L$i"
ncks -v w,rho,theta CAM_6_4_060_06032025.mpasa120_mpasa120.FHS94.mountain_gw.mpasa120_L$i.cam.h0i.0001-01-01-21600.native.nc nCells.nc
echo "Renaming mislabeled dimensions"
ncrename -O -d nCells,ncol nCells.nc
echo "Extracing correctly dimensioned variables"
ncks -x -v w,rho,theta CAM_6_4_060_06032025.mpasa120_mpasa120.FHS94.mountain_gw.mpasa120_L$i.cam.h0i.0001-01-01-21600.native.nc ncol.nc
echo "Concatenating into output file."
ncks -A nCells.nc ncol.nc
ncremap -m ~/regrid/map_mpasa120_to_1x1_aave.nc ncol.nc CAM_6_4_060_06032025.mpasa120_mpasa120.FHS94.mountain_gw.mpasa120_L$i.cam.h0i.0001-01-01-21600.nc
rm -v nCells.nc ncol.nc

done