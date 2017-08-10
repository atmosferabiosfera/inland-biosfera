#!/usr/bin/env python

# =============================================================================
# imports

import sys, os
import numpy

# import netcdf package, either netCDF4 (recommended, support netcdf-4 files) or scipy.io.netcdf
# actually... scipy.io.netcdf is not supported yet...
try:
    from netCDF4 import Dataset, num2date
    use_netcdf4 = True
    use_compression = True
except ImportError:
    use_netcdf4 = False
    use_compression = False
    print('please install netcdf4-python from http://netcdf4-python.googlecode.com')
    sys.exit( 1 )
    #use_netcdf4 = False
    #try:
    #    from scipy.io import netcdf 
    #except ImportError:
    #    print('please install netCDF4 or scipy.io.netcdf')
    #    sys.exit( 1 )
    #try:
    #    from netcdftime import num2date #get this from netCDF4-python if using scipy
    #except ImportError:
    #    print('please install num2date.py from netCDF4')
    #    sys.exit( 1 )

# =============================================================================
# Functions

def Usage():
    print('Usage: inland-retile.py ifile ofile <istile={0,1}> ')
    print('')
    sys.exit( 1 )

# get nodata value for a given variable
def nodataval(var):
    if "missing_value" in var.ncattrs():
        nodata = var.getncattr("missing_value") 
    elif "_FillValue" in var.ncattrs():
        nodata = var.getncattr("_FillValue") 
    else:
        nodata = None
    return nodata

# =============================================================================
# Arguments

if len(sys.argv) < 3:
    Usage()
      
ifile=sys.argv[1]
ofile=sys.argv[2]
if len(sys.argv) == 4:
    istile=sys.argv[3]
else:
    istile=0

if not os.path.exists(ifile):
    print('ifile '+ifile+' not found')
    Usage()

print ('ifile: '+str(ifile))
print ('ofile: '+str(ofile))

nvegtypes=18
vars_ignore = ['date', 'vegtype']

if os.environ.get('INLAND_ICOMPRESSOUT') is not None:
    if os.environ.get('INLAND_ICOMPRESSOUT') == '0':
        use_compression = False

# ==============================================================================
# create output file by copying input file dims, global attributes and vars

# if using netcdf4 (unless INLAND_ICOMPRESSOUT=0 is set)
# define file as netcdf4-classic and use zlib compression in vars 
if use_netcdf4:
    ids = Dataset(ifile, 'r')
    if use_compression:
        ods = Dataset(ofile, 'w', format='NETCDF4_CLASSIC')
    else:
        ods = Dataset(ofile, 'w', format='NETCDF3_CLASSIC')
else:
    ids = netcdf.netcdf_file(ifile, 'r')
    ods = netcdf.netcdf_file(ofile, 'w')

# create dimensions - could first define and then write for more efficiency
dimnames=[]
for dimname in ids.dimensions: 

    #print( dimname )
    dim = ids.dimensions[dimname]

    # rename tile dimension to vegtype and rezise to nvegtypes+1
    dname = dimname
    dlen = len(dim)
    if dimname == "tile":
        dname = "vegtype"
        dlen = nvegtypes+1

    d = ods.createDimension( dname, dlen )
    if dimname in ids.variables:
        var = ids.variables[dimname]
        dims = var.dimensions
        dvar = var[:]            
        if dimname == "tile":
            dvar = range(1,dlen+1)
            dims = "vegtype",
        v = ods.createVariable( dname, var.dtype, dims )
        for a in var.ncattrs():
            v.setncattr( a, var.getncattr(a) )       
        v[:] = dvar
    d = None

#copy global attributes
for a in ids.ncattrs():
    ods.setncattr( a, getattr(ids,a) )

# create variables
for varname in ids.variables:

    if (varname in ids.dimensions) or (varname in vars_ignore):
        continue
    var = ids.variables[varname]

    # rename tile dimension to vegtype
    #dim2 = var.dimensions
    dim2 = []
    for d in list(var.dimensions):
        if d != "tile":
            dim2.append(d)
        else:
            dim2.append('vegtype')
    dim2=tuple(dim2)

    # compress var
    if use_netcdf4 and use_compression:
        v = ods.createVariable( varname, var.dtype, dim2, zlib=True, complevel=2 )
    else:
        v = ods.createVariable( varname, var.dtype, dim2 )

    #add attrs
    for a in var.ncattrs():
        v.setncattr( a, getattr(var,a) )

    v = None


# ==============================================================================
# populate output file by vegtype

# get vegtype map
vegtype_data = ids.variables['vegtype']
numtiles = len(ids.variables["tile"][:])

# process each var, except dims and any other special vars

vars = ods.variables
dimvars = ods.dimensions
# close dataset, so memory doesn't explode when writing vars
ods.close()

for varname in vars:
#for varname in ["npptot"]:

    if (varname in dimvars) or (varname in vars_ignore):
        continue
    print("==== "+varname)

    ivar = ids.variables[varname]
    numdims = len(ivar.shape)

    ods = Dataset(ofile, 'r+')
    ovar = ods.variables[varname]
    # fill with nodata
    nodata = nodataval(ovar)
    if istile == "1":
        nodata = 0
    ovar[:] = nodata
        
    # process each "vegtype"
    for j in range(1,nvegtypes+1):

        data = ovar[:,j-1]

        # process each tile in ifile
        for i in range(0,numtiles-1):
            idata = ivar[:,i]
            # for first tile, inititialize with nodata
            if i == 0:
                data[:] = nodata
                if istile == "1": 
                    data2 = 0
            # get data that exists for vegtype j in tile i
            mask = numpy.equal( vegtype_data[:,i], j )           
            if istile == "1": 
                data2 = data2 + numpy.choose( mask, (0,idata) )
            else:
                data = numpy.choose( mask, (data,idata) )

        if istile == "1": 
            ovar[:,j-1] = data2
        else:
            ovar[:,j-1] = data

    # copy last tile (average)
    ovar[:,nvegtypes] = ivar[:,numtiles-1]

    # close dataset, so memory doesn't explode when writing vars
    ods.close()
    
# ==============================================================================
# close files

ids.close()
#ods.close()


