#!/bin/bash

# this script is to prepare the crops map file for use in inland-agro model

# define input & output parameters here
ofile="crops_monfreda_5min.nc"
ofile2="crops_monfreda_05d.nc"
# directory where 5min crop files are stored
# get them from http://www.geog.mcgill.ca/landuse/pub/Data/175crops2000/NetCDF/
idir="../../../data/monfreada/www.geog.mcgill.ca/landuse/pub/Data/175crops2000/NetCDF/"
# which crops should be included in output file
crops="soybean maize wheat sugarcane rice"
    
rm -f tmp*.nc $ofile ${ofile}.gz $ofile2

# create grid description files
cat > tmp_grid_5min.txt <<EOF
gridtype  = lonlat
xname     = longitude
yname     = latitude
xsize     = 4320
ysize     = 2160
xfirst    = -179.958
xinc      = 0.0833333
yfirst    = 89.9583
yinc      = -0.0833359
EOF

cat > tmp_grid_05d.txt <<EOF
gridtype  = lonlat
xname     = longitude
yname     = latitude
xsize     = 720
ysize     = 360
xfirst    = -179.75
xinc      = 0.5
yfirst    = 89.75
yinc      = -0.5
EOF

# loop for each crop: extract and set level
codes=""
crops=( $crops )
for (( i = 0 ; i < ${#crops[*]} ; i++ )) ; do
    crop=${crops[$i]}
    codes=$codes$(($i+1))" : "${crop}" / "
    ifile=${crop}_5min.nc
    # unzip zip file, don't overwrite
    unzip -n ${idir}/${ifile}.zip
    # extract first level of each files and set output levels
    cdo -O setlevel,$(($i+1)) -sellevel,1 ${ifile} tmp${i}.nc
done

# merge all files to single 5min file
cdo -O merge tmp*.nc tmp_.nc 
cdo -O setmisstoc,0 -setgrid,tmp_grid_5min.txt tmp_.nc $ofile
rm -f tmp*.nc

# interpolate to 05d file
cdo interpolate,tmp_grid_05d.txt $ofile $ofile2

# fix metadata
for f in $ofile $ofile2; do
    ncatted -h -O -a units,level,m,c,"none" $f
    ncatted -h -O -a codes,level,c,c,"${codes}" $f
    ncatted -h -O -a units,cropdata,m,c,"fraction" $f
    ncatted -h -O -a standard_name,longitude,m,c,"longitude" $f
    ncatted -h -O -a units,longitude,m,c,"degrees_east" $f
    ncatted -h -O -a standard_name,latitude,m,c,"latitude" $f
    ncatted -h -O -a units,latitude,m,c,"degrees_north" $f
    ncatted -h -O -a history,global,d,, $f
done


cdo sellonlatbox,-82.0,-34.25,13.1,-55.0 crops_monfreda_05d.nc crops_monfreda_05d_sam.nc

gzip $ofile

rm -f tmp_grid_5min.txt tmp_grid_05d.txt

# plot
#../../inland-contrib/contrib/inland-compare-plot2D.py crops_monfreda_05d.nc img 1 0 1 
#../../inland-contrib/contrib/inland-compare-plot2D.py crops_monfreda_05d_sam.nc img 1 0 1
