#!/usr/bin/env python

import sys, os

import pylab as pl
from array import array
from matplotlib.ticker import * #MultipleLocator, FormatStrFormatter

#from netCDF4 import Dataset, num2date
from scipy.io import netcdf 
from netcdftime import num2date #get this from netCDF4-python


# =============================================================================
# Functions

def Usage():
    print('Usage: plot-1D.py ifile <odir> <plot_2d=0> <plot_norm=0>')
    print('')
    sys.exit( 1 )

def plot_var(varname,ofile,title,units,xaxis,yaxis):
    #print('plot_var('+varname+','+ofile+','+title+',...)')
    pl.clf()
    fig, ax = pl.subplots()
    if plot_norm:
        title = title + ' / norm. diff ( % ) '
    elif units != 'none' and units != 'None':
        title = title + ' ( ' + units + ' )'
    pl.suptitle(title)
    pl.axhline(y=0, color='0.5')
    #pl.plot(xaxis,yaxis,'g-')
    #pl.plot(xaxis,yaxis,'bo')
    pl.plot(xaxis,yaxis,color='gray', linestyle='solid', linewidth=2)
    pl.plot(xaxis,yaxis,'ko')
    xborder=float(xaxis[1]-xaxis[0])/10
    pl.xlim(xaxis[0]-xborder,xaxis[len(xaxis)-1]+xborder)
    numlocs = len(xaxis)-1
    if numlocs > 10:
        numlocs = 7
    # no idea why this check needs to be done...
    if len(xaxis) == 2:
        ax.xaxis.set_major_locator(MultipleLocator( numlocs ) )
    else:
        ax.xaxis.set_major_locator(MaxNLocator( numlocs ))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
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
#enable this to plot 2d (e.g. pft) vars
plot_norm = False
if len(sys.argv)>=5:
    if (sys.argv[4]=='1'):
        plot_norm = True

exclude_vars=[ 'time_weights', 'longitude', 'latitude' ]

print('plot-1D.py '+ifile+' '+ifile_base+' '+odir)

if not os.path.exists(ifile):
    print('ifile '+ifile+' not found')
    Usage()


# =============================================================================
# Main

# open file
#ncfile = Dataset(ifile)
ncfile = netcdf.netcdf_file(ifile, 'r')
times = ncfile.variables['time']
dates = num2date(times[:],units=times.units,calendar=times.calendar)
dates_axis = []#array('i')
for date in dates:
    dates_axis.append(int(str(date)[0:4]))

# loop all vars
for var_name in ncfile.variables: 
    var = ncfile.variables[var_name]
    ndims = len(var.shape)
    ofile_base = os.path.splitext(ifile_base)[0]+'_'+var_name
    if ndims < 3:
        continue
    if var_name in exclude_vars:
        continue
    # plot 1d var (file has 3 dims)
    if ndims == 3:
        ofile_name = odir + '/' + ofile_base + '.png'
        ofile_title = var_name
        plot_var(var_name,ofile_name,ofile_title,var.units,dates_axis,var[:,0,0])
    # plot 2d var (pft)
    elif plot_2d:
        odir2 = odir + '/' + var.dimensions[1]
        if not os.path.exists(odir2):
            os.mkdir(odir2)
        # loop for all z-dims
        for i in range(0,var.shape[1]):
            ofile_name = odir2 + '/' + ofile_base + '_' + str(i).zfill(2) + '.png'
            ofile_title = var_name + ' / ' + var.dimensions[1] + ' = ' + str(i)
            plot_var(var_name,ofile_name,ofile_title,var.units,dates_axis,var[:,i,0,0])

# close file
ncfile.close()
ncfile = None

print('done')

