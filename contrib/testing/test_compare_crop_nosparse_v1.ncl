;Compares output crop variable (sparse) with reference
begin

reffname            =   "ref_output/crop-yearly-2006.nc"
inpfname            =   "../../output/crop-yearly-2006.nc"
vname               =   "cropyld"
indpft              =   13

plistfname          =   "points.list"

arqref = addfile(reffname,"r")
varref = arqref->$vname$

arqinp = addfile(inpfname,"r")
varinp = arqinp->$vname$

printVarSummary(varref)
printVarSummary(varinp)

;Read the points list, using wc to count lines beforehand
npoi = toint(systemfunc("wc -l " + plistfname + "|grep -o '[0-9][0-9]*'"))
print("Number of points (check first and last lines in input file): " + npoi)
plist =  asciiread(plistfname, (/npoi,3/), "float")
pids = plist(:,0)
plons = plist(:,1)
plats = plist(:,2)


print(pids)
print(lats)

npoi = dimsizes(pids)
print("PID       Reference       Input")
do i = 0,npoi-1
    valinp = varinp(0,{indpft},{lats(i)},{lons(i)})
    valref = varref(0,{indpft},{lats(i)},{lons(i)})
    ;print(valinp)
    ;print(valref)
    print(pids(i)+"        "+valref+"      "+valinp)
end do


end