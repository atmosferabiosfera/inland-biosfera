load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;Default values for example and testing purposes, should be commented in order to be set as command arguments in shell
;fname       = "prec.mon.nc"
;vname       = "prec"
;infolder    = "../../input/"
;outfolder   = "../../../clim_input_sparse/"

plistfname  = "points.list" ;lon lat, no headers, no trailing newlines

timvname = "time"
levvname = "level" ;This can change below
latvname = "latitude"
lonvname = "longitude"

;Handle the cases where level is actually layer
vcaseslayer=(/"sandpct","claypct"/)
if (any(vcaseslayer.eq.vname)) then
	levvname="layer"
end if

fin = addfile(infolder + fname,"r")
varin = fin->$vname$
;Also extract coordinate variables for ease of use later
intime = varin&$timvname$
inlevs = varin&$levvname$
inlats = varin&$latvname$
inlons = varin&$lonvname$

printVarSummary(varin)

;Read the points list, using wc to count lines beforehand
npoi = toint(systemfunc("wc -l " + plistfname + "|grep -o '[0-9][0-9]*'"))
print("Number of points (check first and last lines in input file): " + npoi)
plist =  asciiread(plistfname, (/npoi,3/), "float")
pids = plist(:,0)
plons = plist(:,1)
plats = plist(:,2)


;Now,the default structure is that latitude is a coordinate variable, but its just the point ids (should not repeat), and longitude is singleton as well as level when applicable. Actual lat and lon values will go on the auxlatitude and auxlongitude variables, that should also be written in all outputs
;DEPRECATED:The default structure is only latitude is a coordinate variable, and as such it must be sorted, and longitude is singleton as well as level when applicable
;In line with the grid default, latitude is sorted from positive to negative
permvec = dim_pqsort(pids, 1)
pids = pids(permvec)
plats = plats(permvec)
plons = plons(permvec)

;Allocate output variable
varout = new((/dimsizes(intime),dimsizes(inlevs),dimsizes(plats),1/),typeof(varin))
varout!0 = timvname
varout!1 = levvname
varout!2 = latvname
varout!3 = lonvname

;Copy attributes (units, long_name, missing value, etc.)
copy_VarAtts(varin,varout)

;Also copy metadata to coordvars
outtime = intime
outlevs = inlevs
outlats = todouble(pids) ;Point ID's as double
outlons = plons(0) ;Set as just the first point
copy_VarAtts(intime, outtime)
copy_VarAtts(inlevs, outlevs)
copy_VarAtts(inlats, outlats)
copy_VarAtts(inlons, outlons)

;Add a comment in the coordvars' metadata
outlats@sparse_comment = "Sparse array, pids as latitude and singleton dummy longitude. Actual point coordinates are in auxlongitude and auxlatitude"
outlons@sparse_comment = "Sparse array, pids as latitude and singleton dummy longitude. Actual point coordinates are in auxlongitude and auxlatitude"

;Generate auxiliary latitude and longitude variables, not sorted but corresponding to the actual values
auxlats = plats
copy_VarAtts(inlats, auxlats)
auxlons = plons
copy_VarAtts(inlons, auxlons)

;Assign coordinate variables. Remember that longitude is now singleton
varout&$timvname$ = outtime
varout&$levvname$ = outlevs
varout&$latvname$ = outlats
varout&$lonvname$ = outlons

printVarSummary(varout)

do i = 0,npoi-1
    varout(:,:,i,0) = (/ varin(:,:,{plats(i)},{plons(i)}) /) ;Mind the singleton level here
end do


;Erase output file if exists, open it and write variables
system("[ -e" + outfolder + fname + " ] && rm " + outfolder + fname)
fout = addfile(outfolder + fname,"c")

;Copy global attributes and set a new one
copy_VarAtts(fin,fout)
fout@sparse_comment = "Sparse array, pids as latitude and singleton dummy longitude. Actual longitudes are in auxlongitude. Extracted with script " + get_script_prefix_name() + " on " + systemfunc("date")


;Write output and auxiliary variables
fout->$vname$ = varout
fout->auxlatitude = auxlats
fout->auxlongitude = auxlons
fout->auxpid = pids


end
