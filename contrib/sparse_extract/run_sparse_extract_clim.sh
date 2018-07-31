#!/bin/bash
# gabriel.abrahao@ufv.br 
# Runs the NCL script that extracts sparse matrices from gridded data for a whole folder
# Guesses the variable name from the file name using a regex that may have to be changed depending on the application

nclbin="ncl"

nclscript="sparse_extract_clim_v0.2.ncl"

#infolder="../../../clim_input/"
infolder="../../input/"
outfolder="../../../clim_input_sparse/"

#fname="surta.nc"

for fname in $(cd $infolder;ls *.nc); do
	#Guess variable name
	vname=$(echo $fname | sed 's/\([^.]*\)\..*/\1/')

	#Handle clay and sand filenames "soita."
	if [ $vname = "soita" ] ; then
		vname=$(echo $fname | sed 's/\([^.]*\)\.\([^.]*\)\..*/\2/')
	fi


	#Handle clay and sand variable cases where there is an added "pct"
	if [ $vname = "sand" ] || [ $vname = "clay" ] ; then
		vname=$vname"pct"
	fi
	
	echo "============================= RUNNING FILE: "$fname" VARIABLE: "$vname" ================================================="
	$nclbin -Q 'fname="'$fname'"' 'vname="'$vname'"' 'infolder="'$infolder'"' 'outfolder="'$outfolder'"' $nclscript
	
done
