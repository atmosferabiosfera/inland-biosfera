load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

fname = "surta.nc"
vname = "surta"
infolder = "../../../clim_input"
outfolder = "../../../clim_input_sparse/"

plistfname = "points.list"

latvname = "latitude"
lonvname = "longitude"


fin = addfile(fname,"r")
varin = fin->$vname$

printVarSummary(varin)




end