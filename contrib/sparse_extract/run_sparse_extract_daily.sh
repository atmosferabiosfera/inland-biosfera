#!/bin/bash
# gabriel.abrahao@ufv.br 
# Runs the NCL script that extracts sparse matrices from gridded data for a whole folder
# Guesses the variable name from the file name using a regex that may have to be changed depending on the application

nclbin="ncl"

nclscript="sparse_extract_daily_v0.2.ncl"

#infolder="../../../clim_input/"
infolder="/home/gabriel/agroserv/xavier/preproc/daily_all/"
outfolder="../../../daily_input_sparse/"

#Not the smartest approach here, looping though all files every year and running only the ones of the year. But allows for not defining file names or variable prefixes a priori
syear=2005
eyear=2007

#fname="surta.nc"

for fname in $(cd $infolder;ls *.nc); do
	#Guess variable name
	vname=$(echo $fname | sed 's/\([^.]*\)\..*/\1/')

	for year in $(seq $syear $eyear); do
		if [ $(echo $fname | grep -o '[0-9]\{4\}') == $year ]; then
			#echo $fname'	|	'$vname
			
			echo "============================= RUNNING FILE: "$fname" VARIABLE: "$vname" ================================================="
			$nclbin -Q 'fname="'$fname'"' 'vname="'$vname'"' 'infolder="'$infolder'"' 'outfolder="'$outfolder'"' $nclscript
		fi
	done
done
