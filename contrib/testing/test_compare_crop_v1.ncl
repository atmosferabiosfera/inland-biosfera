;Compares output crop variable (sparse) with reference
begin

reffname            =   "ref_output/crop-yearly-2006.nc"
inpfname            =   "../../output/crop-yearly-2006.nc"
vname               =   "cropyld"
indpft              =   13

arqref = addfile(reffname,"r")
varref = arqref->$vname$

arqinp = addfile(inpfname,"r")
varinp = arqinp->$vname$
lats = arqinp->latitude
;Remember that those have 2 dimensions, the singleton longitude is a must for inland writes too
lons = arqinp->auxlongitude(:,0)
pids = arqinp->auxpid(:,0)

printVarSummary(varref)
printVarSummary(varinp)


print(pids)
print(lats)

npoi = dimsizes(pids)
print("PID       Reference       Input")
do i = 0,npoi-1
    valinp = varinp(0,{indpft},i,0)
    valref = varref(0,{indpft},{lats(i)},{lons(i)})
    ;print(valinp)
    ;print(valref)
    print(pids(i)+"        "+valref+"      "+valinp)
end do


end