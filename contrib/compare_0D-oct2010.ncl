; This script generates comparative graphs for every variable on 0D ibis/inland
; as specified (must match the WRITE() command in the code, the column order is
; not automatically inferred as the column names are no longer printed on the
; model's .csv file.
;
; This script can also be used as a base to convert the .csv file into an
; equivalent .nc file (NetCDF).
;
; Important: This version has the column sets designed to work with Hewlley's
; 0D IBIS version dated oct/2010. There's a new version to the apr/2012 output
; format.
;
; author: Fabricio Zimmerer Murta <fabricio@familiamurta.org> 2011-12-23
; Last change: fzm 2012-06-14

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Data file names
inland_file = "inland.c1.k83.lbamip.csv"
ibis_file = "ibis.c1.k83.lbamip"

; Calendar type for data input. "Standard" (which counts leap years), "365_day",
; "360_day", "gregorian", "proleptic_gregorian", "noleap", "julian".
calendar = "365_day"

; determine time interval to plot (format "YYYY-MM-DD HH")
datestart="2010-01-01 00"
dateend = "2010-01-02 00"

; set a filename prefix to differentiate the output
file_prefix="1day"

; Specify global CSV delimiter
delimiter = ","

; CSV Columns array
columns = (/ "iyear", "imonth", "iday", "istep", "jday", "swnet", "lwnet",    \
             "-fvapa*hvap", "-fsena", "soihfl", "delcanheat", "NaN1",         \
             "-fvapa", "grunof", "NaN2", "gdrain", "gdrain+grunof",           \
             "delsoilmoist", "NaN3", "NaN4", "NaN5", "vegT", "tsoi1", "NaN6", \
             "albedo", "NaN7", "fpar", "soilmoist", "soilmoist2",             \
             "soilmoist3", "soilmoist4", "soilmoist5", "soilmoist6",          \
             "tsoi(1,1)", "tsoi(1,2)", "tsoi(1,3)", "tsoi(1,4)", "tsoi(1,5)", \
             "tsoi(1,6)", "SoilWet", "ginvap", "gtrans", "gsuvap", "NaN8",    \
             "NaN9", "CanopInt", "tgpptot*0.012", "tnpptot*0.012",            \
             "-tneetot*0.012", "autoresp*0.012", "tco2mic*0.012", "NaN10",    \
             "NaN11", "NaN12", "cbioltot", "cbiortot", "cbiowtot", "clitll",  \
             "clitlm", "clitls", "clitrl", "clitrm", "clitrs", "clitwl",      \
             "clitwm", "clitws", "totlivbio", "cbiowtot", "totlai" /)

; These are filename-friendly version of the variables above
columnsffd = (/ "iyear", "imonth", "iday", "istep", "jday", "swnet", "lwnet", \
             "neg-fvapa_times_hvap","neg-fsena", "soihfl", "delcanheat","NaN1",\
             "neg-fvapa", "grunof", "NaN2", "gdrain", "gdrain_plus_grunof",   \
             "delsoilmoist", "NaN3", "NaN4", "NaN5", "vegT", "tsoi1", "NaN6", \
             "albedo", "NaN7", "fpar", "soilmoist", "soilmoist2",             \
             "soilmoist3", "soilmoist4", "soilmoist5", "soilmoist6",          \
             "tsoi_1", "tsoi_2", "tsoi_3", "tsoi_4", "tsoi_5",                \
             "tsoi_6", "SoilWet", "ginvap", "gtrans", "gsuvap", "NaN8",       \
             "NaN9", "CanopInt", "tgpptot", "tnpptot",            \
             "neg-tneetot", "autoresp", "tco2mic", "NaN10",    \
             "NaN11", "NaN12", "cbioltot", "cbiortot", "cbiowtot", "clitll",  \
             "clitlm", "clitls", "clitrl", "clitrm", "clitrs", "clitwl",      \
             "clitwm", "clitws", "totlivbio", "cbiowtot", "totlai" /)

; shortcut commonly used formatted units
wm2=":F0:W:F18:V:F0:m:S:-2" ; W.m^-2
jm2=":F0:J:F18:V:F0:m:S:-2" ; J.m^-2
jm2s=jm2+":N,F18:V:F0:s:S:-1" ; ...s^-1
kgm2s=":F0:Kg:F18:V:F0:m:S:-2:N,F18:V:F0:s:S:-1" ; Kg.m^-2.s^-1
kgh2om2=":F0:Kg:B:H:B:2:N:O:N,F18:V:F0:m:S:-2" ; Kg(H2O).m^-2
kgh2om2s=kgh2om2+":N,F18:V:F0:s:S:-1" ; ...s^-1
kelvin="degrees Kelvin"

; The columns measured are *0.012 which collects only the carbon contribuition
; of mass
kgcm2=":F0:Kg:B:C:N,F18:V:F0:m:S:-2" ; Kg(C).m^-2
kgcm2s=kgcm2+":N,F18:V:F0:s:S:-1" ; ...s^-1

; Column units
colunits = (/ "years", "months", "days", "steps", "days", wm2, wm2,           \
             jm2s ,wm2, wm2, jm2, "NaN1",                                     \
             kgm2s, kgh2om2s, "NaN2", kgh2om2s, kgh2om2s,                     \
             kgh2om2, "NaN3", "NaN4", "NaN5", kelvin, kelvin, "NaN6",         \
             "fraction", "NaN7", wm2, kgh2om2, kgh2om2,                       \
             kgh2om2, kgh2om2, kgh2om2, kgh2om2,                              \
             kelvin, kelvin, kelvin, kelvin, kelvin,                          \
             kelvin, "fraction", kgh2om2s, kgh2om2s, kgh2om2s, "NaN8",       \
             "NaN9", "water storage", kgcm2s, kgcm2s,            \
             kgcm2s, kgcm2s, kgcm2s, "NaN10",    \
             "NaN11", "NaN12", kgcm2, kgcm2, kgcm2, kgcm2,  \
             kgcm2, kgcm2, kgcm2, kgcm2, kgcm2, kgcm2,      \
             kgcm2, kgcm2, kgcm2, kgcm2, "fraction" /)

; set the option to go with cd_inv_calendar into getsteps, getstartsteps
calopt = 0
calopt@calendar = calendar

; Get default missing value from a double. We will update this later.
defmisval = new(1, double)
missing_value = defmisval@_FillValue
delete(defmisval)

function get_trim_str(handle, column)
begin
 return(stringtoint(str_strip(str_get_field(handle,column,delimiter))))
end

function get_data_fld(handle, column)
local i,data
begin
 data = stringtodouble(str_get_field(handle,column,delimiter))
 data@_FillValue = missing_value
; Turn NaN's into missing values.
 if (any(isnan_ieee(data))) then
  replace_ieeenan(data,data@_FillValue,0)
 end if
 return(data)
end

function indexof(value)
local index,result
begin
 result = -1
 do index = 0, dimsizes(columns)-1
  if (columns(index) .eq. value) then
   result = index
  end if
 end do
 return result
end

; Get number of steps to run based on ISO YYYY-MM-DD HH start and end.
function getsteps()
local interval, year, month, day, hour, nil, steps
begin
 interval = (/ datestart, dateend /)
 year = stringtoint(str_get_field(interval,1,"-"))
 month = stringtoint(str_get_field(interval,2,"-"))
 day = stringtoint(str_get_field(interval,3,"-"))
 hour = stringtoint(str_get_field(interval,2," "))
 nil = hour-hour
 steps = cd_inv_calendar(year,month,day,hour,nil,nil, \
                         "hours since "+year(0)+"-1-1 00:00:00",calopt)
 return(toint(steps(1)-steps(0)))
end

function getstartstep(units)
local year, month, day, hour, nil
begin
 year = stringtoint(str_get_field(datestart,1,"-"))
 month = stringtoint(str_get_field(datestart,2,"-"))
 day = stringtoint(str_get_field(datestart,3,"-"))
 hour = stringtoint(str_get_field(datestart,2," "))
 nil = hour-hour ; this has to be set here to guarantee it is at ymdh dimsize
 return(toint(cd_inv_calendar(year,month,day,hour,nil,nil,units,calopt)))
end

function isNaN(column)
local haystack
begin
 haystack=str_split_by_length(columns(column),3)
 if (haystack(0).eq."NaN") then
  return(True)
 else
  return(False)
 end if
end

begin
 inland = 0
 ibis = 1

 datafile = new( (/2/), string)
 datafile(inland) = inland_file
 datafile(ibis) = ibis_file

; Read data
 print("Reading data.")
 data_inland = asciiread(datafile(inland),-1,"string")
 data_ibis = asciiread(datafile(ibis),-1,"string")
 rawdata = (/ data_inland, data_ibis /)

; Determine missing value for variables based on a NaN column.
 missing_value = get_data_fld(rawdata(0,0),indexof("NaN1")+1)

 print("Preparing plot.")

; Sanity checks
 nrows = dimsizes(rawdata(0,:))
 ncols = dimsizes(str_split(rawdata(0,0),delimiter))
 if (ncols .ne. dimsizes(str_split(rawdata(1,0), delimiter))) then
  print("ERROR: column count mismatches between compared files!")
  status_exit(1)
 end if
 if (ncols .ne. dimsizes(columns)) then
  print("ERROR: column count mismatches between known columns and data files!")
  status_exit(1)
 end if

 year = get_trim_str(rawdata,1)
 month = get_trim_str(rawdata,2)
 day = get_trim_str(rawdata,3)
 hour = get_trim_str(rawdata,4)-1 ; hour in datafile is 1-24
  
 nil = hour-hour
 timestep = cd_inv_calendar(year,month,day,hour,nil,nil, \
                            "hours since "+min(year(:,0))+"-1-1 00:00:00",0)

; Make up a nice date label
 timecharlocal = sprinti("%0.2i",day)+"/"+sprinti("%0.2i",month)+"/"+ year + \
                 ":C:"+sprinti("%0.2i",hour)+"::00::00 UTC"
; timecharlocal = day + "/" + month + "/" + year + ":C:" + hour + "::00::00 UTC"

 res = True
 res@gsnScale = True 
 res@gsnMaximize = True
 res@vpHeightF = 0.43 ;0.43
 res@vpWidthF = 0.85  ;0.85

 res@gsnDraw = False
 res@gsnFrame = False
 res@tmLabelAutoStride = True ; do not let tick labels overlap

; Determine Panel properties
 panel_res = True
 panel_res@gsnScale = True
 panel_res@gsnMaximize = True

; Deduce initial and final steps
 steps = getsteps()
 fstep = getstartstep(timestep@units) ; first step
 lstep = fstep+steps                  ; last step

; Sanity checks on start and end steps
 time_adjusted = False
 if (steps.gt.dimsizes(rawdata(0,:))) then
  print("WARNING: Number of time steps requested is greater than available time steps on data. Reverting to run the full data interval.")
  steps = dimsizes(rawdata(0,:))
  fstep = 0
  lstep = steps-1
  time_adjusted = True
 end if
 if (fstep.lt.0) then
  print("WARNING: first step was set to before data is available. Fixing.")
  steps = steps + fstep ; something + (-something) to keep end date
  fstep = 0
  time_adjusted = True
 end if
 if (lstep.ge.dimsizes(rawdata(0,:))) then
  print("WARNING: Last step set beyond last data is available. Fixing.")
  lstep = dimsizes(rawdata(0,:))-1
  steps = lstep - fstep
  time_adjusted = True
 end if

 ;res@tmXBFormat = "f" ; remove trailing .0 from date labels
 res@trXMinF = timestep(0,fstep)
 res@trXMaxF = timestep(0,lstep)

; Map X indexes from ugly timesteps to nice datetime strings. :)
 res@tmXBMode = "Explicit"
 res@tmXBValues = timestep(0,fstep:lstep:steps/4) ;pick up 5 values
 res@tmXBLabels = timecharlocal(0,fstep:lstep:steps/4)
 res@tmXBMinorOn = True
 res@tmXBMinorValues = timestep(0,fstep:lstep:steps/20) ;pick up 5 values
 res@tiXAxisString = ":F0:time"
 print("Interval "+timecharlocal(0,fstep)+" to "+timecharlocal(0,lstep))
 if (time_adjusted) then
  print("*** NOTICE: Time interval adjusted to fit the available data!")
 end if

; Define workstation default options
 worksttn_props = "png"
 worksttn_props@wkWidth = "800"
 worksttn_props@wkHeight = "800"

; NOTE: If you want to make a transparent background image
; 1. Write a .ps file instead of a .png (make the page in landscape)
; 2. Convert it to png with imagemagick's convert tool:
;    convert -transparent white file.ps file.png
 plot = new(2, graphic)
 do column=indexof("swnet"), dimsizes(columns)-1 - indexof("swnet")

; Process only non-NaN columns. :)
  if (.not. isNaN(column)) then
   worksttn_props@wkFileName=sprinti("%0.2i-",column)
   if (.not. str_is_blank(file_prefix)) then
    worksttn_props@wkFileName=worksttn_props@wkFileName+file_prefix+"-"
   end if
   worksttn_props@wkFileName=worksttn_props@wkFileName+columnsffd(column)+".png"
   print("Graphs col#"+column+": "+columns(column)+" / "+ worksttn_props@wkFileName)

; Get Y axis units label from units array
   res@tiYAxisString = colunits(column)

; Show legends on this first graph
   res@pmLegendDisplayMode = "Always"
   res@pmLegendHeightF = 0.10
; Adjust legend names with data names
   res@xyExplicitLegendLabels = new(2, string)
   res@xyExplicitLegendLabels(ibis) = "IBIS 0D/Hewlley - oct/2010"
   res@xyExplicitLegendLabels(inland) = "InLand 0D/rev68f352a9f600"


; Define lines color and thickness for first graph
   res@xyLineThicknesses = (/2.0, 2.0/)
   res@xyLineColors = (/ "blue", "green" /)

; Choose title, filename, set workstation and PLOT!
   res@tiMainString = columns(column)+":: 0D InLand and IBIS"
   worksttn = gsn_open_wks(worksttn_props, "0D InLand and IBIS")
   plot(0) = gsn_csm_xy(worksttn, timestep(:,fstep:lstep),              \
                        get_data_fld(rawdata(:,fstep:lstep),column+1),res)

; Prevent legends from being shown on this graph (but create object)
   res@pmLegendDisplayMode = "Never"

; Reset lines color and thickness for the second graph
   delete(res@xyLineThicknesses)
   delete(res@xyLineColors)

; Choose title, filename, refresh workstation and PLOT!
   res@tiMainString = columns(column)+":: 0D IBIS - InLand difference"

   plot(1) = gsn_csm_xy(worksttn, timestep(:,fstep:lstep),                     \
                        get_data_fld(rawdata(ibis,fstep:lstep),column+1) -     \
                        get_data_fld(rawdata(inland,fstep:lstep),column+1),res)

; Make a fancy panel containing both graphs
   gsn_panel(worksttn, plot, (/2,1/), panel_res)
  end if
 end do
end
