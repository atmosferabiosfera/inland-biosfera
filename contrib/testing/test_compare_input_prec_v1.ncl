begin
;Compares output crop variable (sparse) with reference

reffname            =   "../../input/daily/prec.daily.2005.nc"
inpfname            =   "../../output/inland-daily-2005.nc"
refvname            =   "prec"
inpvname            =   "rain"


arqref = addfile(reffname,"r")
varref = arqref->$refvname$

arqinp = addfile(inpfname,"r")
varinp = arqinp->$inpvname$
lats = arqinp->latitude
;Remember that those have 2 dimensions, the singleton longitude is a must for inland writes too
lons = arqinp->auxlongitude(:,0)
pids = arqinp->auxpid(:,0)

printVarSummary(varref)
printVarSummary(varinp)


print(pids)
print(lats)

npoi = dimsizes(pids)
print("PID      Day       Reference       Input     Diff.")
do i = 0,npoi-1
    do day = 0,dimsizes(varref&time)-1
        valinp = varinp(day,i,0)
        valref = varref(day,{lats(i)},{lons(i)})
        diff = valinp - valref
        ;print(valinp)
        ;print(valref)
        print(pids(i)+"     "+day+"        "+valref+"      "+valinp+"       "+diff)
    end do
end do


end