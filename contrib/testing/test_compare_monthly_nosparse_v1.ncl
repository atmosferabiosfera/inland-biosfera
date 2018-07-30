;Compares output crop variable (sparse) with reference
begin

reffname            =   "ref_output/inland-monthly-2006.nc"
inpfname            =   "../../output/inland-monthly-2006.nc"
vname               =   "lail"
;indpft              =   13
indmon              =   1

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
print(plats)

npoi = dimsizes(pids)
print("MONTH = " + indmon)
print("PID       Reference       Input")
do i = 0,npoi-1
    ;valinp = varinp(indmon,{indpft},{plats(i)},{(i)})
    ;valref = varref(indmon,{indpft},{plats(i)},{(i)})
    valinp = varinp(indmon,{plats(i)},{plons(i)})
    valref = varref(indmon,{plats(i)},{plons(i)})
    ;print(valinp)
    ;print(valref)
    print(pids(i)+"        "+valref+"      "+valinp)
end do


end