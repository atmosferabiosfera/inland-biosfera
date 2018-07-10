;Compares output crop variable (sparse) with reference
begin

reffname            =   "ref_output/crop-yearly-2006.nc"
inpfname            =   "../../output/crop-yearly-2006.nc"
vname               =   "cropyld"


arqref = addfile(reffname,"r")
varref = arqref->$vname$

arqinp = addfile(inpfname,"r")
varinp = arqinp->$vname$
lats = arqinp->latitude
lons = arqinp->auxlongitude
pids = arqinp->auxpid

printVarSummary(varref)
printVarSummary(varinp)

end