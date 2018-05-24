#!/bin/bash

usage="
usage:   inland-maketiles.sh ifile ofile #tiles <tilefrac_levels> <cdo_commands> 
         where #tiles+1 (last tile is water) values each for tilefrac_levels and cdo_commands (use "" for none)
         if last cdo_command is not specified then '-setrtoc,0,200,30' is used

example: inland-maketiles.sh vegtype.nc bla1.nc 2 0.65 0.25 0.10 '' 'setvals,1,9,11,14' ''
"

:<<COMMENT

This script is used to generate vegtype files with tiles

example usage:

# these commands create files with a second tile where forest(1) and cerrado(9) are replaced by pasture(18) in varying proportions (10, 25 and 50%) 
inland-maketiles.sh vegtype.nc vegtype-deforest-10p.nc 2 0.90 0.10 0.0 "" "setvals,1,18,9,18" ""
inland-maketiles.sh vegtype.nc vegtype-deforest-25p.nc 2 0.75 0.25 0.0 "" "setvals,1,18,9,18" ""
inland-maketiles.sh vegtype.nc vegtype-deforest-50p.nc 2 0.50 0.50 0.0 "" "setvals,1,18,9,18" ""

# a simple example with 2 identical tiles and second tile has 0 tilefrac values (used as a template)
inland-maketiles.sh vegtype.nc vegtype-mlpt02-equal.nc 2 1 0 0 "" "" ""
# another simple example with 2 different tiles and 10% waterfrac
inland-maketiles.sh vegtype.nc vegtype-mlpt02-diff.nc 2 0.75 0.25 0 "" "setvals,1,9,9,11,11,14" ""

The same can be acheived using these commands

# a simple example with 2 identical tiles and second tile has 0 tilefrac values (used as a template)

cd tiles
mkdir tmp
rm -f tmp/vegtype-t0?.nc tmp/tilefrac-t0?.nc
# create vegtype levels and modify tiles
cdo setlevel,1 ../vegtype.nc tmp/vegtype-t01.nc
cdo setlevel,2 ../vegtype.nc tmp/vegtype-t02.nc
cdo  setlevel,3 -setrtoc,0,200,30 ../vegtype.nc tmp/vegtype-t03.nc
# modify any of the tilefrac levels, here is a simple example
cdo setlevel,1 -setrtoc,0,100,1 -chvar,vegtype,tilefrac ../vegtype.nc tmp/tilefrac-t01.nc 
cdo setlevel,2 -setrtoc,0,100,0 -chvar,vegtype,tilefrac ../vegtype.nc tmp/tilefrac-t02.nc 
cdo setlevel,3 -setrtoc,0,100,0 -chvar,vegtype,tilefrac ../vegtype.nc tmp/tilefrac-t03.nc 
# merge files
cdo -O merge tmp/vegtype-t0?.nc tmp/vegtype-mlpt02.nc
cdo -O merge tmp/tilefrac-t0?.nc tmp/tilefrac-mlpt02.nc
cdo -O merge tmp/vegtype-mlpt02.nc tmp/tilefrac-mlpt02.nc vegtype-mlpt02-equal.nc
cd ..

# another simple example with 2 different tiles and 10% waterfrac

cd tiles
mkdir tmp
rm -f tmp/vegtype-t0?.nc tmp/tilefrac-t0?.nc
# create vegtype levels and modify tiles, here is a simple example where forest->savanna, savanna->shrubs and shrubs->desert
cdo setlevel,1 ../vegtype.nc tmp/vegtype-t01.nc
cdo setlevel,2 -setvals,1,9,9,11,11,14 ../vegtype.nc tmp/vegtype-t02.nc
cdo  setlevel,3 -setrtoc,0,200,30 ../vegtype.nc tmp/vegtype-t03.nc
# modify any of the tilefrac levels, here is a simple example
cdo setlevel,1 -setrtoc,0,100,0.75 -chvar,vegtype,tilefrac ../vegtype.nc tmp/tilefrac-t01.nc 
cdo setlevel,2 -setrtoc,0,100,0.25 -chvar,vegtype,tilefrac ../vegtype.nc tmp/tilefrac-t02.nc 
cdo setlevel,3 -setrtoc,0,100,0 -chvar,vegtype,tilefrac ../vegtype.nc tmp/tilefrac-t03.nc 
# merge files
cdo -O merge tmp/vegtype-t0?.nc tmp/vegtype-mlpt02.nc
cdo -O merge tmp/tilefrac-t0?.nc tmp/tilefrac-mlpt02.nc
cdo -O merge tmp/vegtype-mlpt02.nc tmp/tilefrac-mlpt02.nc vegtype-mlpt02-diff.nc
cd ..

COMMENT

#set -x
#echo $# - "$*"

# make sure min. args.
if [[ $# -lt 3 ]] ; then echo "$usage" ; exit ; fi

# parse min args
ifile=$1
ifile_base=`basename $1 .nc`
ofile=$2
numtiles=$3

if [[ ! -f $ifile ]] ; then echo -e "\nifile $ifile missing..." ; echo "$usage" ; exit ; fi
if [[ $numtiles -lt 2 ]] ; then echo -e "\n#tiles must be > 1 " ; echo "$usage" ; exit ; fi

# make sure args count ok
target=$((3+($numtiles+1)*2))
if [[ $# -lt $(($target-1)) ]] ; then echo "$usage" ; exit ; fi
if [[ $# -gt $target ]] ; then echo "$usage" ; exit ; fi
#if [[ $# -lt $target ]] ; then args=("${@}" "") ; else args=("${@}") ; fi
args=("${@}")
if [[ ${args[$(($target-1))]} == "" ]] ; then args[$(($target-1))]="setrtoc,0,200,30" ; fi
#echo ${#args[@]} - ${args[@]}

cd tiles
if [[ ! -d tmp ]] ; then mkdir tmp ; fi
rm -f tmp/vegtype-t??.nc tmp/tilefrac-t??.nc

totfrac=0

for (( i = 1 ; i <= (($numtiles+1)) ; i++ )) ; do

ii=`printf %02d $i`
#echo $i - $ii - ${args[$((3+($numtiles+1)+$i-1))]}

# create vegtype levels and modify tiles, here is a simple example where forest->savanna and shrubs->desert
cdo_command=${args[$((3+($numtiles)+$i))]}
if [[ $cdo_command != "" ]] ; then cdo_command="-"$cdo_command ; fi
cdo -s setlevel,${i} $cdo_command ../$ifile tmp/vegtype-t${ii}.nc

# modify any of the tilefrac levels, here is a simple example
frac=${args[$((2+$i))]}
totfrac=`echo $totfrac+$frac | bc`
if [[ $frac == "" ]] ; then frac=0 ; fi

cdo -s setlevel,${i} -setrtoc,0,100,${frac} -chvar,vegtype,tilefrac ../${ifile} tmp/tilefrac-t${ii}.nc

done

# make sure total fraction count = 1, if not warn user
if [ $totfrac != "1" -a $totfrac != "1.0" -a $totfrac != "1.00" ] ; then echo -e "\nWARNING: fraction total = $totfrac , expecting 1.0\n" ; fi

# merge files
ii=`printf %02d $i`
cdo -s -O merge tmp/vegtype-t??.nc tmp/vegtype-txx.nc
cdo -s -O merge tmp/tilefrac-t??.nc tmp/tilefrac-txx.nc
cdo -s -O merge tmp/vegtype-txx.nc tmp/tilefrac-txx.nc $ofile

# cleanup
rm -f tmp/vegtype-t??.nc tmp/tilefrac-t??.nc


cd ..
