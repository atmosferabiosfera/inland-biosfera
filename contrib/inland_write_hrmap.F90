#include "inland_config.h"

module inland_write_hrmap

contains


! this subroutine creates a HR map for each variables in ifile, in function of the
! HR map hrfile and tile map tfile
subroutine whrmap (ifile, hrfile, tfile, ofile)

  use inland_netcdfutils ! from src/ dir

  implicit none

! -----------------------------------------------------------------------
! input variables
  character*(*) :: ifile, hrfile, tfile, ofile

! -----------------------------------------------------------------------
! local variables
  integer ids, hrds, tds, ods
  integer i, j, k, p, c, ierr, natts, j1, j2
  integer ndims, dimid, dimlen, dimidtile
  integer nvars, varid, varlen, vartype
  integer dimids(5)
  integer istart(5),icount(5)
  integer chunksize(5)
  character*80 dimname, varname, attname
  character*1024 tmpchar
  integer ntiles, ntilemap, npfts, tmpnpfts
  integer nodataparent, nodatachild
  real*8 nodata
  integer icompressout

  real*8, dimension(:), allocatable :: dimdata
  real*8, dimension(:,:,:), allocatable :: idata
  real*8, dimension(:), allocatable :: odata
  integer, dimension(:,:,:), allocatable :: itilechild
  integer, dimension(:), allocatable :: ihrtileparent
  integer, dimension(:,:), allocatable :: tilemap
  integer nlon1, nlat1, nlon2, nlat2, ndimdata

#ifdef GFORTRAN
#include <netcdf.inc>
#else /* GFORTRAN */
  include 'netcdf.inc'
#endif

  data istart / 1,1,1,1,1 /, icount / 1,1,1,1,1 /

  dimidtile = -1
  npfts = 1
  ndimdata = 1
  
  print *,ifile, ' ', hrfile, ' ', tfile, ' ', ofile

  ! open datasets for reading
  print *,"opening datasets"
  ierr = NF_OPEN(ifile, NF_NOWRITE, ids)
  if (ierr .ne. NF_NOERR) call netcdf_error("error opening "//ifile,ierr)
  ierr = NF_OPEN(hrfile, NF_NOWRITE, hrds)
  if (ierr .ne. NF_NOERR) call netcdf_error("error opening "//hrfile,ierr)
  ierr = NF_OPEN(tfile, NF_NOWRITE, tds)
  if (ierr .ne. NF_NOERR) call netcdf_error("error opening "//tfile,ierr)

  ! create output dataset
  ! must use netcdf-4 because of lon/lat size exceeds netcdf-classic limits
  ! also compress to save on file size
  icompressout = 1
  if ( icompressout .eq. 0 ) then
     ierr = NF_CREATE(ofile,NF_CLOBBER,ods)
  else
#ifdef NETCDF4
     ierr = NF_CREATE(ofile,OR(NF_CLOBBER, OR(nf_netcdf4, nf_classic_model)),ods)
#else
     ierr = NF_CREATE(ofile,NF_CLOBBER,ods)
     print *, 'Warning in whrmap netcdf-4 not available so you cannot use compression'
#endif
  end if
  if (ierr .ne. NF_NOERR) call netcdf_error("error opening "//ofile,4)
  
  !print *, ids, ' ', hrds, ' ', tds, ' ', ods

  ! define dimensions
  ierr = NF_INQ_NDIMS(ids, ndims)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting ndims "//ifile,ierr)
  print *,"creating",ndims,"dims"
  do i = 1, ndims
     ierr = NF_INQ_DIM(ids, i, dimname, dimlen)
     if (ierr .ne. NF_NOERR) call netcdf_error("error getting dimname "//ifile,ierr)
     !resize lon and lat to hrmap
     if ( ( dimname .eq. "longitude" ) .or. ( dimname .eq. "latitude") ) then
        ierr = NF_INQ_DIMID(hrds, dimname, dimid)
        if (ierr .ne. NF_NOERR) call netcdf_error("error getting dimid "//hrfile,ierr)
        ierr = NF_INQ_DIMLEN(hrds, dimid, dimlen)
        if (ierr .ne. NF_NOERR) call netcdf_error("error getting dimlen "//hrfile,ierr)
        print *,i, trim(dimname), dimlen
     end if
     ndimdata = MAX(ndimdata, dimlen)
     print *,i, trim(dimname), dimlen
     if ( dimname .eq. "tile" ) dimidtile = i
     if ( dimname .eq. "pft" ) npfts = dimlen
     ierr = NF_DEF_DIM(ods, dimname, dimlen, dimid)
     if (ierr .ne. NF_NOERR) call netcdf_error("error creating dimname "//dimname//" for "//ofile,ierr)    
  end do

  ! define variables
  ierr = NF_INQ_NVARS(ids, nvars)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting nvars for "//ifile,ierr)
  print *,"creating",nvars,"vars"
  do i = 1, nvars
     ! get dimids from ifile - should make sure these match in ofile...
     ierr = NF_INQ_VAR(ids, i, varname, vartype, ndims, dimids, natts)
     if (ierr .ne. NF_NOERR) call netcdf_error("error getting varname "//ifile,ierr)
     ! skip tile dimension
     ! TODO test with pft
     if ( ndims .gt. 2 ) then
        do j = 1, ndims-1
           if ( dimids(j) .ge. dimidtile ) dimids(j) = dimids(j+1)
        end do
        ndims = ndims - 1
     end if
     print *,i, trim(varname), ndims, dimids(1:ndims)
     ! define var
     ierr = NF_DEF_VAR(ods, varname, vartype, ndims, dimids, varid)
     if (ierr .ne. NF_NOERR) call netcdf_error("error creating varname "//trim(varname)//" for "//ofile,ierr)

     if ( icompressout .gt. 0 ) then
        ierr = NF_DEF_VAR_DEFLATE(ods,varid,0,1,icompressout)
        if (ierr .ne. NF_NOERR) then
           print *, 'Warning in inivar at NF_DEF_VAR_DEFLATE ',varname
           print *, NF_STRERROR(ierr)
        end if
      end if

     ! add attributes
     call netcdf_copyattrs(varname,natts,ids,i,ods,varid)
  end do !nvars

  ! add global attributes 
  ierr = NF_INQ_NATTS(ids, natts)
  call netcdf_copyattrs('NF_GLOBAL',natts,ids,NF_GLOBAL,ods,NF_GLOBAL)

  ! exit define mode
  ierr = NF_ENDDEF(ods)
  if (ierr .ne. NF_NOERR) call netcdf_error("error enddef for"//ofile,ierr)

  ! copy dimension vars from ifile to ofile, except lon/lat
  allocate(dimdata(ndimdata))
  ierr = NF_INQ_NDIMS(ids, ndims)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting ndims "//ifile,ierr)
  print *,"copying",ndims,"dims from ifile"
  do i = 1, ndims
     ierr = NF_INQ_DIM(ids, i, dimname, dimlen)
     if (ierr .ne. NF_NOERR) call netcdf_error("error getting dim "//ifile,ierr)

     ! skip lon/lat
     if ( (dimname .eq. "longitude" ) .or.(dimname .eq. "latitude" ) ) cycle

     ! read data
     ierr = NF_INQ_VARID(ids, dimname, varid)
     if (ierr .ne. NF_NOERR) call netcdf_error("error getting varid "//ifile,ierr)
     istart(:) = 1 
     icount(:) = 1
     icount(1) = dimlen
     ierr = NF_GET_VARA_DOUBLE(ids, varid, istart, icount, dimdata)
     if (ierr .ne. NF_NOERR) call netcdf_error("error getting "//trim(dimname)//" data in "//ifile,ierr) 

     ! write data
     ierr = NF_INQ_VARID(ods, dimname, varid)
     if (ierr .ne. NF_NOERR) call netcdf_error("error getting varid "//ofile,ierr)
     ierr = NF_PUT_VARA_DOUBLE(ods, varid, istart, icount, dimdata)
     if (ierr .ne. NF_NOERR) call netcdf_error("error putting "//trim(varname)//" data in "//ofile,ierr)
  end do

  ! copy lon/lat dimension vars from hrfile to ofile - duplicate code from above...
  ierr = NF_INQ_NDIMS(hrds, ndims)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting ndims "//hrfile,ierr)
  print *,"copying lon/lat dims from hrfile"
  do i = 1, ndims
     ierr = NF_INQ_DIM(hrds, i, dimname, dimlen)
     if (ierr .ne. NF_NOERR) call netcdf_error("error getting dim "//hrfile,ierr)

     ! skip other vars
     if ( (dimname .ne. "longitude" ) .and. (dimname .ne. "latitude" ) ) cycle

     ! read data
     ierr = NF_INQ_VARID(hrds, dimname, varid)
     if (ierr .ne. NF_NOERR) call netcdf_error("error getting varid "//hrfile,ierr)
     istart(:) = 1 
     icount(:) = 1
     icount(1) = dimlen
     ierr = NF_GET_VARA_DOUBLE(hrds, varid, istart, icount, dimdata)
     if (ierr .ne. NF_NOERR) call netcdf_error("error getting "//trim(dimname)//" data in "//hrfile,ierr) 

     ! write data
     ierr = NF_INQ_VARID(ods, dimname, varid)
     if (ierr .ne. NF_NOERR) call netcdf_error("error getting varid "//ofile,ierr)
     ierr = NF_PUT_VARA_DOUBLE(ods, varid, istart, icount, dimdata)
     if (ierr .ne. NF_NOERR) call netcdf_error("error putting "//trim(varname)//" data in "//ofile,ierr)
  end do
  

  ! get tile map

  ! first get lon/lat of ifile and hrfile
  ierr = NF_INQ_DIMID(ids, 'longitude', dimid)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting dimid "//ifile,ierr)
  ierr = NF_INQ_DIMLEN(ids, dimid, nlon1)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting dimlen "//ifile,ierr)
  ierr = NF_INQ_DIMID(ids, 'latitude', dimid)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting dimid "//ifile,ierr)
  ierr = NF_INQ_DIMLEN(ids, dimid, nlat1)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting dimlen "//ifile,ierr)

  ierr = NF_INQ_DIMID(hrds, 'longitude', dimid)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting dimid "//hrfile,ierr)
  ierr = NF_INQ_DIMLEN(hrds, dimid, nlon2)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting dimlen "//hrfile,ierr)
  ierr = NF_INQ_DIMID(hrds, 'latitude', dimid)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting dimid "//hrfile,ierr)
  ierr = NF_INQ_DIMLEN(hrds, dimid, nlat2)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting dimlen "//hrfile,ierr)

  varname = 'tile'
  ierr = NF_INQ_DIMID(tds, varname, dimid)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting "//trim(varname)//" dimid "//tfile,ierr)
  ierr = NF_INQ_DIMLEN(tds, dimid, ntiles)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting ntiles "//tfile,ierr)
  varname = 'itilechild'
  ierr = NF_INQ_VARID(tds, varname, varid)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting "//trim(varname)//" varid "//tfile,ierr)
  allocate(itilechild(nlon1,nlat1,ntiles))
  icount(1) = nlon1
  icount(2) = nlat1
  icount(3) = ntiles
  ierr = NF_GET_VARA_INT(tds, varid, istart, icount, itilechild)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting "//trim(varname)//" var "//tfile,ierr) 
  nodatachild = get_nodata_int(tds, varid)
  
  allocate(ihrtileparent(nlon2 * nlat2))
  varname = 'ihrtileparent'
  ierr = NF_INQ_VARID(hrds, varname, varid)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting "//trim(varname)//" varid "//tfile,ierr)
  icount(1) = nlon2
  icount(2) = nlat2
  icount(3) = 1
  istart(:) = 1
  ierr = NF_GET_VARA_INT(hrds, varid, istart, icount, ihrtileparent(:))
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting "//trim(varname)//" var "//tfile,ierr)
  nodataparent = get_nodata_int(hrds, varid)

  ntilemap = MAXVAL(itilechild)
  allocate(tilemap(3,ntilemap))
  tilemap(:,:) = nodataparent
  do k = 1, ntiles
     do j = 1, nlat1
        do i = 1, nlon1
           c = itilechild(i,j,k)
           !print *,c,nodatachild,nodataparent
           ! skip nodata points
           if ( c .eq. nodatachild ) cycle
           tilemap(1,c) = i
           tilemap(2,c) = j
           tilemap(3,c) = k
        end do
     end do
  end do

  ! process each var, except dims and any other special vars

  print *,"filling vars"
  allocate(idata(nlon1,nlat1,ntiles))
  allocate(odata(nlon2 * nlat2))

  ierr = NF_INQ_NVARS(ods, nvars)
  if (ierr .ne. NF_NOERR) call netcdf_error("error getting nvars for"//ofile,ierr)
  do i = 1, nvars

     ! get variable info
     ierr = NF_INQ_VAR(ods, i, varname, vartype, ndims, dimids, natts)
     if (ierr .ne. NF_NOERR) call netcdf_error("error getting varname "//ofile,ierr)

     ! skip dimension variables
     ierr = NF_INQ_DIMID(ods, varname, dimid)
     if (ierr .eq. NF_NOERR) cycle

     print *,i,trim(varname),ndims

     tmpnpfts = 1
     if ( ndims .gt. 3 ) tmpnpfts = 12

     ! get nodata for this variable
     if ( ( vartype .eq. NF_FLOAT ) .or. ( vartype .eq. NF_DOUBLE ) ) then
        nodata = get_nodata_real(ods, i)
     else
        nodata = real(get_nodata_int(ods, i))
     end if

     ! loop over pfts (if needed)
     do p = 1, tmpnpfts

        if ( ndims .gt. 3 ) print *,"pft: ",p

        ! fill odata with it nodata
        odata(:) = nodata

        ! read input variable into idata
        icount(1) = nlon1
        icount(2) = nlat1
        if ( ndims .gt. 3 ) then
           istart(3) = p
           icount(3) = 1
           istart(4) = 1
           icount(4) = ntiles
        else
           istart(3) = 1
           icount(3) = ntiles
           istart(4) = 1
           icount(4) = 1
        end if

        ierr = NF_INQ_VARID(ids, varname, varid)
        if (ierr .ne. NF_NOERR) call netcdf_error("error getting "//trim(varname)//" varid "//ifile,ierr)
        ierr = NF_GET_VARA_DOUBLE(ids, varid, istart, icount, idata)
        if (ierr .ne. NF_NOERR) call netcdf_error("error getting "//trim(varname)//" data in "//ifile,ierr)
        
        ! loop over all points in odata
        do j = 1, nlon2 * nlat2
           c = ihrtileparent(j)
           ! if point has no parent tile, skip it
           if ( c .eq. nodataparent ) then
              cycle
           end if
           odata(j) = idata(tilemap(1,c),tilemap(2,c),tilemap(3,c))
        end do
        
        ! write odata to ofile
        icount(1) = nlon2
        icount(2) = nlat2
        icount(3) = 1
        istart(4) = 1
        icount(4) = 1
        if ( ndims .gt. 3 ) then
           istart(3) = p
        else
           istart(3) = 1
        end if
        ierr = NF_INQ_VARID(ods, varname, varid)
        if (ierr .ne. NF_NOERR) call netcdf_error("error getting "//trim(varname)//" varid "//ofile,ierr)
        ierr = NF_PUT_VARA_DOUBLE(ods, i, istart, icount, odata)
        if (ierr .ne. NF_NOERR) call netcdf_error("error putting "//trim(varname)//" data in "//ofile,ierr)
        
     end do ! pfts
     
  end do ! nvars

  ! deallocate memory

  deallocate(dimdata)
  deallocate(itilechild)
  deallocate(ihrtileparent)
  deallocate(tilemap)
  deallocate(idata)
  deallocate(odata)

  ! close files
  ierr = NF_CLOSE(ids)
  if (ierr .ne. NF_NOERR) call netcdf_error("error closing "//ifile,ierr)    
  ierr = NF_CLOSE(hrds)
  if (ierr .ne. NF_NOERR) call netcdf_error("error closing "//hrfile,ierr)    
  ierr = NF_CLOSE(tds)
  if (ierr .ne. NF_NOERR) call netcdf_error("error closing "//tfile,ierr)    
  ierr = NF_CLOSE(ods)
  if (ierr .ne. NF_NOERR) call netcdf_error("error closing "//ofile,ierr)    

end subroutine whrmap

end module inland_write_hrmap

! main program
program main

  use inland_write_hrmap

  implicit none

  integer i, iargc, numarg
  character*255 odir, ifile, hrfile, tfile, ofile 

  ! process command-line args
  numarg = iargc()
  if ( numarg .ne. 4 ) then
     print *,""
     print *,"Usage: inland-write--hrmap ifile.nc hrfile.nc tfile.nc ofile.nc"
     print *,""
     stop
  end if

  call getarg ( 1, ifile )
  call getarg ( 2, hrfile )
  call getarg ( 3, tfile )
  call getarg ( 4, ofile )

  ! run
  call whrmap( ifile, hrfile, tfile, ofile)

  print *,'done'

  stop

end program main

