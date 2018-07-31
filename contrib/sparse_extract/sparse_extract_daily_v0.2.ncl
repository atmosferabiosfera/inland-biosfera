load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;Currently orks only with rddaydims=1, i.e. daily data in 3d (lat,lon,time)

;Default values for example and testing purposes, should be commented in order to be set as command arguments in shell
; fname       = "temp.daily.2005.nc"
; vname       = "temp"
; infolder    = "/home/gabriel/agroserv/xavier/preproc/daily_all/"
; outfolder   = "../../../daily_input_sparse/"

plistfname  = "points.list" ;pid lon lat, no headers, no trailing newline but last line must have a newline at its end (check with vim :set list)

timvname = "time"
levvname = "level" ;TODO support 4D variables for use with rddaydims=0
latvname = "latitude"
lonvname = "longitude"


fin = addfile(infolder + fname,"r")
varin = fin->$vname$
;Also extract coordinate variables for ease of use later
intime = varin&$timvname$
;inlevs = varin&$levvname$
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


;The default structure is only latitude is a coordinate variable, and as such it must be sorted, and longitude is singleton as well as level when applicable
;In line with the grid default, latitude is sorted from positive to negative
permvec = dim_pqsort(plats, -1)
pids = pids(permvec)
plats = plats(permvec)
plons = plons(permvec)

print(pids)
print(plats)
print(plons)

;Allocate output variable
;varout = new((/dimsizes(intime),dimsizes(inlevs),dimsizes(plats),1/),typeof(varin)) ;TODO support 4D variables for use with rddaydims=0
varout = new((/dimsizes(intime),dimsizes(plats),1/),typeof(varin))
varout!0 = timvname
;varout!1 = levvname
varout!1 = latvname
varout!2 = lonvname

;Copy attributes (units, long_name, missing value, etc.)
copy_VarAtts(varin,varout)

;Also copy metadata to coordvars
outtime = intime
;outlevs = inlevs
;outlats = todouble(pids) ;Point ID's as double
outlats = plats
outlons = plons(0) ;Set as just the first point
copy_VarAtts(intime, outtime)
;copy_VarAtts(inlevs, outlevs)
copy_VarAtts(inlats, outlats)
copy_VarAtts(inlons, outlons)

;Add a comment in the coordvars' metadata
outlats@sparse_comment = "Sparse array, sorted by latitude and singleton dummy longitude. Actual point coordinates are in auxlongitude and auxlatitude"
outlons@sparse_comment = "Sparse array, sorted by latitude and singleton dummy longitude. Actual point coordinates are in auxlongitude and auxlatitude"

;Generate auxiliary latitude, longitude and pid variables, not sorted but corresponding to the actual values
auxlats = plats
copy_VarAtts(inlats, auxlats)
auxlons = plons
copy_VarAtts(inlons, auxlons)
auxpids = pids



;Assign coordinate variables. Remember that longitude is now singleton
varout&$timvname$ = outtime
;varout&$levvname$ = outlevs
varout&$latvname$ = outlats
varout&$lonvname$ = outlons

printVarSummary(varout)

do i = 0,npoi-1
    varout(:,i,0) = (/ varin(:,{plats(i)},{plons(i)}) /) ;Mind the absence of a singleton level here ;TODO support 4D variables for use with rddaydims=0
end do


;Erase output file if exists, open it and write variables
system("[ -e" + outfolder + fname + " ] && rm " + outfolder + fname)
fout = addfile(outfolder + fname,"c")

;Copy global attributes and set a new one
copy_VarAtts(fin,fout)
fout@sparse_comment = "Sparse array, sorted by latitude and singleton dummy longitude. Actual longitudes are in auxlongitude. Extracted with script " + get_script_prefix_name() + " on " + systemfunc("date")


;Write output and auxiliary variables
fout->$vname$ = varout
fout->auxlatitude = auxlats
fout->auxlongitude = auxlons
fout->auxpid = pids


end
