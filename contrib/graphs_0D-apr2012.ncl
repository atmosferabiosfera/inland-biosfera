; This script generates comparative graphs for every variable on 0D ibis/inland
; as specified (must match the WRITE() command in the code, the column order is
; not automatically inferred as the column names are no longer printed on the
; model's .csv file.
;
; This script can also be used as a base to convert the .csv file into an
; equivalent .nc file (NetCDF).
;
; author: Fabricio Zimmerer Murta <fabricio@familiamurta.org> 2011-12-23
; Last change: fzm 2012-04-20

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Data file names
inland_file = "single_point-output.csv"

; Number of heading lines (skip these lines for headers, for example)
headlines = 1

; Calendar type for data input. "Standard" (which counts leap years), "365_day",
; "360_day", "gregorian", "proleptic_gregorian", "noleap", "julian".
calendar = "365_day"

; the output file has only the timesteps on the first column. Here you can
; specify the first year reference so the date can be inferred from the hourly
; timestep on the data.
refyear = 2001

; determine time interval to plot (format "YYYY-MM-DD HH")
datestart="2001-01-01 00"
dateend = "2001-01-01 23"

; set a filename prefix to differentiate the output
file_prefix="1day"

; Specify global CSV delimiter
delimiter = " "

; CSV Columns array
columns = (/ "step", "swnet", "lwnet", "-fvapa*hvap", "-fsena", "soihfl",     \
             "delcanheat", "NaN1", "-fvapa", "grunof", "NaN2", "gdrain",      \
	     "gdrain+grunof", "delsoilmoist", "NaN3", "NaN4", "NaN5", "vegT", \
	     "tsoi1", "NaN6", "albedo", "NaN7", "fpar", "soilmoist1",         \
	     "soilmoist2", "soilmoist3", "soilmoist4", "soilmoist5",          \
	     "soilmoist6", "tsoi(1,1)", "tsoi(1,2)", "tsoi(1,3)", "tsoi(1,4)",\
	     "tsoi(1,5)", "tsoi(1,6)", "SoilWet", "ginvap", "gtrans",         \
	     "gsuvap", "NaN8", "NaN9", "CanopInt", "tgpptot*0.012",           \
	     "tnpptot*0.012", "-tneetot*0.012", "autoresp*0.012",             \
	     "tco2mic*0.012", "ayanpptot", "NaN10", "falll", "cbioltot",      \
	     "cbiortot", "cbiowtot", "clitll", "clitlm", "clitls", "clitrl",  \
	     "clitrm", "clitrs", "clitwl", "clitwm", "clitws", "agb",         \
	     "totlivbio", "cbiowtot", "totlai", "xta+273.15", "xqa", "xua",   \
	     "xprec/1800", "NaN11", "xsin", "xlin", "NaN12" /)

; These are filename-friendly version of the variables above
columnsffd = (/ "step", "swnet", "lwnet", "qle", "qh", "qg", "delcanheat",    \
                "NaN1", "evap", "qs", "NaN2", "qsb", "qt", "delsoilmoist",    \
		"NaN3", "NaN4", "NaN5", "vegt", "baresoilt", "NaN6", "albedo",\
		"NaN7", "fpar", "soilmoist1", "soilmoist2", "soilmoist3",     \
		"soilmoist4", "soilmoist5", "soilmoist6", "tsoi1", "tsoi2",   \
		"tsoi3", "tsoi4", "tsoi5", "tsoi6", "soilwet", "ecanop",      \
		"tveg", "esoil", "NaN8", "NaN9", "CanopInt", "gpp", "npp",    \
		"nee", "autoresp", "heteroresp", "anpp", "NaN10",             \
		"litterfall", "carbpools1", "carbpools2", "carbpools3",       \
		"carbpools4", "carbpools5", "carbpools6", "carbpools7",       \
		"carbpools8", "carbpools9", "carbpools10", "carbpools11",     \
		"carbpools12", "carbpools13", "totlivbiom", "AbvGrndWood",    \
		"lai", "tair", "qair", "wind", "rainf", "NaN11", "swdown",    \
		"lwdown", "NaN12" /)

; shortcut commonly used formatted units
wm2=":F0:W:F18:V:F0:m:S:-2" ; W.m^-2
jm2=":F0:J:F18:V:F0:m:S:-2" ; J.m^-2
jm2s=jm2+":N,F18:V:F0:s:S:-1" ; ...s^-1
kgm2s=":F0:Kg:F18:V:F0:m:S:-2:N,F18:V:F0:s:S:-1" ; Kg.m^-2.s^-1
kgh2om2=":F0:Kg:B:H:B:2:N:O:N,F18:V:F0:m:S:-2" ; Kg(H2O).m^-2
kgh2om2s=kgh2om2+":N,F18:V:F0:s:S:-1" ; ...s^-1
kgh2oair=":F0:Kg:B:H:B:2:N:O:N,F18:V:F0:Kg:B:air:N,S:-1"; Kg(H2O).Kg(air)^-1
kelvin="degrees Kelvin"
mps=":F0:m:F18:V:F0:s:S:-1" ; m.s^-1

; The columns measured are *0.012 which collects only the carbon contribuition
; of mass
kgcm2=":F0:Kg:B:C:N,F18:V:F0:m:S:-2" ; Kg(C).m^-2
kgcm2s=kgcm2+":N,F18:V:F0:s:S:-1" ; ...s^-1

; Column units
colunits = (/ "hours", wm2, wm2, jm2s ,wm2, wm2, jm2, "NaN1", kgm2s, kgh2om2s,\
              "NaN2", kgh2om2s, kgh2om2s, kgh2om2, "NaN3", "NaN4", "NaN5",    \
	      kelvin, kelvin, "NaN6", "fraction", "NaN7", "fraction", kgh2om2,\
	      kgh2om2, kgh2om2, kgh2om2, kgh2om2, kgh2om2, kelvin, kelvin,    \
	      kelvin, kelvin, kelvin, kelvin, "fraction", kgh2om2s, kgh2om2s, \
	      kgh2om2s, "NaN8", "NaN9", "units", kgcm2s, kgcm2s, kgcm2s,      \
	      kgcm2s, kgcm2s, kgcm2s, "NaN10", kgcm2, kgcm2, kgcm2, kgcm2,    \
	      kgcm2, kgcm2, kgcm2, kgcm2, kgcm2, kgcm2, kgcm2, kgcm2, kgcm2,  \
	      kgcm2, kgcm2, kgcm2, "fraction", kelvin, kgh2oair, mps, kgh2om2s,\
	      "NaN11", wm2, wm2 /)

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

; Get number of steps (lines on file)  to run based on ISO YYYY-MM-DD HH start 
; and end.
function getsteps()
;local interval, year, month, day, hour, nil, steps
local interval, steps
begin
 interval = (/ datestart, dateend /)
 year = stringtoint(str_get_field(interval,1,"-"))  ; !!-..-.. ..
 month = stringtoint(str_get_field(interval,2,"-")) ; ..-!!-.. ..
 day = stringtoint(str_get_field(interval,3,"-"))   ; ..-..-!! ..
 hour = stringtoint(str_get_field(interval,2," "))  ; ..-..-.. !!
 nil = hour-hour ; I do not remember why I chose to do this way. :)
 steps = cd_inv_calendar(year,month,day,hour,nil,nil, \
                         "hours since "+year(0)+"-1-1 00:00:00",calopt)
 return(toint(steps(1)-steps(0)))
end

; gets first step to be processed (first line on file)
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
 datafile = inland_file

; Read data
 print("Reading data.")
 data = asciiread(datafile,-1,"string")

; Make rawdata with inland+ibis data skipping the heading lines
 rawdata = data(headlines::)

; If we have just one column, it must mean the delimiter is wrong.
 if (str_fields_count(rawdata(0), delimiter) .le. 1)  then
  print("ERROR: Detected just one column on data files. Maybe the field delimiter is wrong?")
  status_exit(1)
 end if
  
; Determine missing value for variables based on a NaN column.
 missing_value = get_data_fld(rawdata(0),indexof("NaN1")+1)
 print("Determined missing value from NaN1 column as "+missing_value)

 print("Preparing plot.")

; Sanity checks
 ;nrows = dimsizes(rawdata(0,:))
 ncols = str_fields_count(rawdata(0),delimiter)
 if (ncols .ne. dimsizes(columns)) then
  print("ERROR: column count mismatches between known columns and data files!")
  status_exit(1)
 end if

 ; now the reading is inversed. Instead of having year-jday-hour on the data
 ; file, we have just the timestep. Then from the timestep and base start year,
 ; we create the year-month-day-hour information.
 timestep = toint(get_trim_str(rawdata,1))
 timestep@units = "hours since "+refyear+"-01-01 00:00:00"
 utc_date = cd_calendar(timestep, -5)
 year = utc_date(:,0)
 month = utc_date(:,1)
 day = utc_date(:,2)
 hour = utc_date(:,3)
  
 nil = hour-hour

 ; when the data file had year-jday-hour information, we read these informations
 ; directly from the file and turned into timestep with the command below, we
 ; did not need to have a 'reference year' (refyear) to base with, as the first
 ; line of the data file contained the first year of the run.
 ;timestep = cd_inv_calendar(year,month,day,hour,nil,nil, \
 ;                           "hours since "+min(year(:,0))+"-1-1 00:00:00",0)

; Make up a nice date label
 timecharlocal = sprinti("%0.2i",day)+"/"+sprinti("%0.2i",month)+"/"+ year + \
                 ":C:"+sprinti("%0.2i",hour)+"::00::00 UTC"
; timecharlocal = day + "/" + month + "/" + year + ":C:" + hour + "::00::00 UTC"

 res = True
 res@gsnScale = True 
 res@gsnMaximize = True
 res@vpHeightF = 0.43 ;0.43
 res@vpWidthF = 0.85  ;0.85

 res@gsnDraw = True
 res@gsnFrame = True
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
 if (steps.gt.dimsizes(rawdata)) then
  print("WARNING: Number of time steps requested is greater than available time steps on data. Reverting to run the full data interval.")
  steps = dimsizes(rawdata)
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
 if (lstep.ge.dimsizes(rawdata)) then
  print("WARNING: Last step set beyond last data is available. Fixing.")
  lstep = dimsizes(rawdata)-1
  steps = lstep - fstep
  time_adjusted = True
 end if

 ;res@tmXBFormat = "f" ; remove trailing .0 from date labels
 res@trXMinF = timestep(fstep)
 res@trXMaxF = timestep(lstep)

; Map X indexes from ugly timesteps to nice datetime strings. :)
 res@tmXBMode = "Explicit"
 res@tmXBValues = timestep(fstep:lstep:steps/4) ;pick up 5 values
 res@tmXBLabels = timecharlocal(fstep:lstep:steps/4)
 res@tmXBMinorOn = True
 res@tmXBMinorValues = timestep(fstep:lstep:steps/20) ;pick up 5 values
 res@tiXAxisString = ":F0:time"
 print("Interval: "+timecharlocal(fstep)+" to "+timecharlocal(lstep))
 if (time_adjusted) then
  print("*** NOTICE: Time interval adjusted to fit the available data!")
 end if

; Define workstation default options
 worksttn_props = "png"
 worksttn_props@wkWidth = "800"
 worksttn_props@wkHeight = "600"

; NOTE: If you want to make a transparent background image
; 1. Write a .ps file instead of a .png (make the page in landscape)
; 2. Convert it to png with imagemagick's convert tool:
;    convert -transparent white file.ps file.png
; plot = new(1, graphic)
 do column=indexof("swnet"), dimsizes(columns)-1 - indexof("swnet")

; Process only non-NaN columns. :)
  if (.not. isNaN(column)) then
   worksttn_props@wkFileName=sprinti("%0.2i-",column)
   if (.not. str_is_blank(file_prefix)) then
    worksttn_props@wkFileName=worksttn_props@wkFileName+file_prefix+"-"
   end if
   worksttn_props@wkFileName=worksttn_props@wkFileName+columnsffd(column)+".png"
   print("Graph col#"+column+": "+columns(column)+" / "+ worksttn_props@wkFileName)

; Get Y axis units label from units array
   res@tiYAxisString = colunits(column)

; Show legends on this first graph
   res@pmLegendDisplayMode = "Always"
   res@pmLegendHeightF = 0.10
; Adjust legend names with data names
   res@xyExplicitLegendLabels = "InLand 0D"

; Define lines color and thickness for first graph
   res@xyLineThicknesses = (/2.0/)
   res@xyLineColors = (/ "blue" /)

; Choose title, filename, set workstation and PLOT!
   res@tiMainString = columns(column)+":: 0D InLand"
   worksttn = gsn_open_wks(worksttn_props, "0D InLand")
   plot = gsn_csm_xy(worksttn, timestep(fstep:lstep),              \
                     get_data_fld(rawdata(fstep:lstep),column+1),res)
  end if
 end do
end
