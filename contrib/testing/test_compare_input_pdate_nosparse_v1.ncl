;Compares output crop variable (sparse) with reference
begin

reffname            =   "../../input/crop_params/plantdoy.2005.nc"
inpfname            =   "../../output/crop-yearly-2006.nc"
refvname            =   "plantdoy"
inpvname            =   "plantdate"
indpft              =   13

plistfname          =   "points.list"

arqref = addfile(reffname,"r")
varref = arqref->$refvname$

arqinp = addfile(inpfname,"r")
varinp = arqinp->$inpvname$

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
print("PID       Reference       Input")
do i = 0,npoi-1
    valinp = varinp(0,{indpft},{plats(i)},{plons(i)})
    valref = varref({plats(i)},{plons(i)})
    ;print(valinp)
    ;print(valref)
    print(pids(i)+"        "+valref+"      "+valinp)
end do


end