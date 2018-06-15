load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

fname       = "surta.nc"
vname       = "surta"
infolder    = "../../../clim_input/"
outfolder   = "../../../clim_input_sparse/"

plistfname  = "points.list" ;lon lat, no headers, no trailing newlines

timvname = "time"
latvname = "latitude"
lonvname = "longitude"

fin = addfile(infolder + fname,"r")
varin = fin->$vname$
;Also extract coordinate variables for ease of use later
intime = varin&$timvname$
inlats = varin&$latvname$
inlons = varin&$lonvname$

printVarSummary(varin)

;Read the points list, using wc to count lines beforehand
npoi = toint(systemfunc("wc -l " + plistfname + "|grep -o '[0-9][0-9]*'"))
print("Number of points (check first and last lines in input file): " + npoi)
plist =  asciiread(plistfname, (/npoi,2/), "float")
plons = plist(:,0)
plats = plist(:,1)


;The default structure is only latitude is a coordinate variable, and as such it must be sorted, and longitude is singleton
;In line with the grid default, latitude is sorted from positive to negative
permvec = dim_pqsort(plats, -1)
plats = plats(permvec)
plons = plons(permvec)

;Allocate output variable
varout = new((/dimsizes(intime),dimsizes(plats),1/),typeof(varin))
varout!0 = timvname
varout!1 = latvname
varout!2 = lonvname

;Also copy metadata to coordvars
outtime = intime
outlats = plats
outlons = plons(0) ;Set as just the first point
copy_VarAtts(intime, outtime)
copy_VarAtts(inlats, outlats)
copy_VarAtts(inlons, outlons)

;Generate an auxiliary longitude variable, not sorted but corresponding to the actual values
auxlons = plons
copy_VarAtts(inlons, auxlons)

;Assign coordinate variables. Remember that longitude is now singleton
varout&$timvname$ = outtime
varout&$latvname$ = outlats
varout&$lonvname$ = outlons

printVarSummary(varout)

do i = 0,npoi-1
    varout(:,i,0) = (/ varin(:,0,{plats(i)},{plons(i)}) /) ;Mind the singleton level here
end do


;Erase output file if exists, open it and write variables
system("[ -e" + outfolder + fname + " ] && rm " + outfolder + fname)
fout = addfile(outfolder + fname,"c")
fout->$vname$ = varout
fout->auxlongitude = auxlons


end