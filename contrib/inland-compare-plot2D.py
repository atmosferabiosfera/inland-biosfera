#!/usr/bin/env python

import sys, os

import pylab as pl
import numpy as np
from mpl_toolkits.basemap import Basemap

use_scipy=True
#use_scipy=False
if use_scipy:
    from scipy.io import netcdf 
    from netcdftime import num2date #get this from netCDF4-python
else:
    from netCDF4 import Dataset, num2date, date2num

# =============================================================================
# Functions

def Usage():
    print('Usage: plot-2D.py ifile <odir> <plot_2d=0> <plot_norm=0>')
    print('')
    sys.exit( 1 )

def plot_var(varname,ofile,title,units,data,fill_value,lon,lat):
    pl.clf()
    if plot_norm:
        title = title + ' / norm. diff ( % ) '
        units = '%'
    elif units == 'none' or units == 'None':
        units = ''
    else:
        title = title + ' ( ' + units + ' )'
        units = '\n'+units
    pl.suptitle(title)
    # fix stupid scipy.io.netcdf _FillValue bug
    if use_scipy:
        data2 = np.copy( data )
        data2[data2==fill_value] = None
    else:
        data2 = data

    lon2=[ lon[0], lon[len(lon)-1] ]
    lat2=[ lat[0], lat[len(lat)-1] ] 
    lonstep=10
    latstep=10
    labelstyle=None
    fontsize=12
    if lon2 == [-179.75, 179.75] and lat2 == [89.75, -89.75]:
        lat2 = [85,-60]
        lonstep=30
        latstep=20
        #labelstyle="+/-"
        fontsize=9

    # TODO fix top/bottom margins when plotting global or sam map...
    map = Basemap(llcrnrlon=min(lon2),llcrnrlat=min(lat2),urcrnrlon=max(lon2),urcrnrlat=max(lat2),projection='mill')
    map.drawcoastlines(linewidth=1.25)
    map.drawparallels(np.arange(round(min(lat2)),max(lat2),latstep),labels=[1,0,0,0],labelstyle=labelstyle,fontsize=fontsize)
    map.drawmeridians(np.arange(round(min(lon2)),max(lon2),lonstep),labels=[0,0,0,1],labelstyle=labelstyle,fontsize=fontsize)
    data3 = map.transform_scalar(np.flipud(data2), lon, np.flipud(lat), len(lon), len(lat))
    # show data
    #cs = pl.contourf(data)
    cs = map.imshow(data3,interpolation='nearest')
    cbar = map.colorbar(cs)
    cbar.ax.set_xlabel(units,ha='left')   

    pl.savefig(ofile)


# =============================================================================
# Arguments

if len(sys.argv) < 2:
    print('got '+str(len(sys.argv)))
    Usage()
      
ifile=sys.argv[1]  
(odir,ifile_base) = os.path.split(ifile)
if len(sys.argv)>=3:
    odir=odir+'/'+sys.argv[2]
plot_2d = False
if len(sys.argv)>=4:
    if (sys.argv[3]=='1'):
        plot_2d = True
plot_norm = False
if len(sys.argv)>=5:
    if (sys.argv[4]=='1'):
        plot_norm = True


plot_4d = False
exclude_vars=[ 'time_weights', 'longitude', 'latitude' ]


print('plot-2D.py '+ifile+' '+ifile_base+' '+odir)

if not os.path.exists(ifile):
    print('ifile '+ifile+' not found')
    Usage()

if not os.path.exists(odir):
    print('odir '+odir+' not found')
    Usage()

# =============================================================================
# Main

# open file
if use_scipy:
    ncfile = netcdf.netcdf_file(ifile, 'r')
else:
    ncfile = Dataset(ifile)
times = ncfile.variables['time']
dates = num2date(times[:],units=times.units,calendar=times.calendar)

firstvar=True
for var_name in sorted(ncfile.variables.iterkeys()): 
    var = ncfile.variables[var_name]
    ndims = len(var.shape)
    ofile_base = os.path.splitext(ifile_base)[0]+'_'+var_name
    if ndims < 3:
        continue    
    if var_name in exclude_vars:
        continue

    # TODO test this with n>5
    step = 1
    if var.shape[0] > 50:
        step = 20
    steps = range(0,var.shape[0],step)
    if not var.shape[0]-1 in steps:
        steps.append(var.shape[0]-1)
    
    if firstvar:
        firstvar=False
        if len(steps) > 0:
            print('steps: ',len(steps),str(steps))
#    print(var_name)
    print var_name,
    sys.stdout.flush()

    for j in steps:
        #print str(j), #not ok w/ python 3!
        if ( len(steps) > 1 ):
            tmp1 = str(dates[j])[0:4] # this only works for yearly files
            tmp2 = '_' + tmp1
            tmp3 = ' - ' + tmp1
        else:
            tmp1 = ''
            tmp2 = ''
            tmp3 = ''
        if ndims == 3:
            ofile_name = odir + '/' + ofile_base + tmp2 + '.png'
            ofile_title = var_name + tmp3
            plot_var(var_name,ofile_name,ofile_title,var.units,var[j],var._FillValue,ncfile.variables['longitude'].data,ncfile.variables['latitude'].data)
        elif plot_2d:
            odir2 = odir + '/' + var.dimensions[1]
            if not os.path.exists(odir2):
                os.mkdir(odir2)
            for i in range(0,var.shape[1]):
                ofile_name = odir2 + '/' + ofile_base + '_' + tmp1 + '_' + str(i).zfill(2) + '.png'
                ofile_title = var_name + ' / ' + var.dimensions[1] + ' = ' + str(i) + ' / ' + tmp1
                plot_var(var_name,ofile_name,ofile_title,var.units,var[j][i],var._FillValue,ncfile.variables['longitude'].data,ncfile.variables['latitude'].data)
print('')

            
ncfile.close()
ncfile = None

print('done')

