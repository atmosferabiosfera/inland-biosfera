;Gabriel Abrahao
;Makes a surta mask based on a base climate file by setting to 1 only pixels that are not missing
begin

clifname		=	"~/agroserv/xavier/preproc/daily_all/prec.daily.2015.nc"
clivname		=	"prec"

reffname		=	"../input/surta.nc"
refvname		=	"surta"

outfname		=	"surta.nc"
outvname		=	"surta"

fcli = addfile(clifname,"r")
varcli = fcli->$clivname$

fref = addfile(reffname,"r")
varref = fref->$refvname$

printVarSummary(varcli)
printVarSummary(varref)

ndimcli = dimsizes(dimsizes(varcli))

if ndimcli.eq.3 then
	varref = where(.not.ismissing(

end
