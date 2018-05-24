#include "inland_config.h"

#ifdef SINGLE_POINT_MODEL
#error "This subroutine should NOT be compiled for 0D INLAND model option."
#endif

! if --enable-highres_map was not passed to configure, do not compile
#ifdef HRMAP

!.....................................................................
subroutine rdhrvegmap(filen,varname)

      use inland_subgrid
      use inland_combcs
      use inland_control, only: env_debug

      implicit none

#ifdef GFORTRAN
#include <netcdf.inc>
#else /* GFORTRAN */
      include 'netcdf.inc'
#endif

! INPUT
! filen - character*(*) - file name from which to read data
! varname - character*(*) - name of variable from which to read
      character*(*) filen, varname

      integer :: nlons,nlats,ierror,tmpint,tmpcount,domtile,scaling
      integer :: i,j,k,inpoi,ilpt,ib,jb,isubgrid
      real*8 :: xres1,yres1,tmpreal
      integer :: idies, idvar, ndims, idlonvar, idlatvar, idlondim, idlatdim, id3d, ierr, itype
      character*80 tmpname
      integer istart(4), icount(4) ! for reading vars
      integer, dimension (1) :: max_loc
      logical ok
      real*8 :: minlon, maxlon, minlat,maxlat

      real*8, dimension(:), allocatable :: lonvalues, latvalues, xinvegbak
      integer, dimension(:,:), allocatable :: values
      integer, dimension(:,:), allocatable :: vegtypecount
      integer, dimension(:), allocatable :: tmpvegtypecount
      integer watervegtype
      real*8 scaling2

      data istart / 1,1,1,1 /, icount / 1,1,1,1 /
      
      ierror = 0

      watervegtype = 30

      print *,'INFO: processing high-res map                 '//trim(filen)

! sanity checks
! -----------------------------------------

!      if ( nvegtype .ne. 15 ) then
!         print *, 'Error in rdhrvegmap, nvegtype != 15, = ',nvegtype
!         print *, NF_STRERROR(ierr)
!         ierror = -1
!         stop 1
!      end if

! Open file
! -----------------------------------------
      ierr = NF_OPEN(filen,NF_NOWRITE,idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in rdhrvegmap, file='//trim(filen)
         print *, NF_STRERROR(ierr)
         ierror = -1
         stop 1
      end if

! Get lon and lat ids
! ------------------
      tmpname='longitude'
      ierr = NF_INQ_VARID(idies,tmpname,idlonvar)
      if (ierr .ne. NF_NOERR) then 
         tmpname='lon'
         ierr = NF_INQ_VARID(idies,tmpname,idlonvar)
         if (ierr .ne. NF_NOERR) then 
            print *, 'Error in rdhrvegmap, cannot get longitude var'
            print *, NF_STRERROR(ierr)
            ierror = -1
            stop 1
         end if
      end if
      ierr = NF_INQ_DIMID(idies,tmpname,idlondim)
      if (ierr .ne. NF_NOERR) then 
         print *, 'Error in rdhrvegmap, cannot get longitude dim'
         print *, NF_STRERROR(ierr)
         ierror = -1
         stop 1
      end if

      tmpname='latitude'
      ierr = NF_INQ_VARID(idies,tmpname,idlatvar)
      if (ierr .ne. NF_NOERR) then 
         tmpname='lat'
         ierr = NF_INQ_VARID(idies,tmpname,idlatvar)
         if (ierr .ne. NF_NOERR) then 
            print *, 'Error in rdhrvegmap, cannot get latitude var'
            print *, NF_STRERROR(ierr)
            ierror = -1
            stop 1
         end if
      end if
      ierr = NF_INQ_DIMID(idies,tmpname,idlatdim)
      if (ierr .ne. NF_NOERR) then 
         print *, 'Error in rdhrvegmap, cannot get latitude dim'
         print *, NF_STRERROR(ierr)
         ierror = -1
         stop 1
      end if


! Get lon and lat length + res (assuming spacing is equal)
! ------------------
      ierr = NF_INQ_DIMLEN(idies,idlondim,nlons)
      if (ierr .ne. NF_NOERR) then 
         print *, 'Error in rdhrvegmap, cannot get longitude length'
         print *, NF_STRERROR(ierr)
         ierror = -1
         stop 1
      end if
      allocate(lonvalues(1:nlons)) 
      icount(1) = nlons
      ierr = NF_GET_VARA_DOUBLE(idies,idlonvar,istart,icount,lonvalues)
      if (ierr .ne. NF_NOERR) then 
         print *, 'Error in rdhrvegmap, cannot get longitude values'
         print *, NF_STRERROR(ierr)
         ierror = -1
         stop 1
      end if
      !print *,lonvalues
      xres1 = lonvalues(2)-lonvalues(1)

      ierr = NF_INQ_DIMLEN(idies,idlatdim,nlats)
      if (ierr .ne. NF_NOERR) then 
         print *, 'Error in rdhrvegmap, cannot get latitude length'
         print *, NF_STRERROR(ierr)
         ierror = -1
         stop 1
      end if 
      allocate(latvalues(1:nlats)) 
      icount(1) = nlats
      ierr = NF_GET_VARA_DOUBLE(idies,idlatvar,istart,icount,latvalues)
      if (ierr .ne. NF_NOERR) then 
         print *, 'Error in rdhrvegmap, cannot get latitude values'
         print *, NF_STRERROR(ierr)
         ierror = -1
         stop 1
      end if
      yres1 = latvalues(2)-latvalues(1)

      if (abs(abs(xres1)-abs(yres1)) .gt. abs(xres1)*10) then 
         print *, 'Error in rdhrvegmap, xres1 != yres1'
         print *, NF_STRERROR(ierr)
         ierror = -1
         stop 1
      end if
      if (abs(abs(xres)-abs(yres)) .gt. abs(xres)*10) then 
         print *, 'Error in rdhrvegmap, xres != yres'
         print *, NF_STRERROR(ierr)
         ierror = -1
         stop 1
      end if


! Make sure HR file encloses domain
! TODO - if it doesn't, fill with global map values instead
! ------------------
      minlon = minval(lonscale(lonindex(1:npoi1))) - xres/2.0
      maxlon = maxval(lonscale(lonindex(1:npoi1))) + xres/2.0
      minlat = minval(latscale(latindex(1:npoi1))) - yres/2.0
      maxlat = maxval(latscale(latindex(1:npoi1))) + yres/2.0
      if ( lonvalues(1) .gt. minlon ) then
         print *, 'Error in rdhrvegmap, W',lonvalues(1),'>',minlon
         ierror = -1
         stop 1
      end if
      if ( lonvalues(nlons) .lt. maxlon ) then
         print *, 'Error in rdhrvegmap, E',lonvalues(nlons),'<',maxlon
         ierror = -1
         stop 1
      end if
      ! latitude may be inverted...
      if ( latvalues(1) .lt. maxlat) then
         print *, 'Error in rdhrvegmap, N',latvalues(1),'<',maxlat
         ierror = -1
         stop 1
      end if
      if ( latvalues(nlats) .gt. minlat ) then
         print *, 'Error in rdhrvegmap, N',latvalues(nlats),'>',minlat
         ierror = -1
         stop 1
      end if

! Calculate portion of HR file to read
! ----------
      tmpreal = minlon
      do i=1,nlons
         if ( lonvalues(i) .gt. tmpreal ) then
            if ( i .eq. 1 ) then
               istart(1) = 1
            else
               istart(1) = i
            end if
            exit
         end if
      end do
      tmpreal = maxlon
      do i=istart(1),nlons
         if ( lonvalues(i) .gt. tmpreal ) then
            icount(1) = i - istart(1)
            exit
         end if
      end do
      if ( icount(1) .eq. 1 ) then
         print *, 'Error in rdhrvegmap, icount(1)=1'
         ierror = -1
         stop 1
      end if
      
      tmpreal = maxlat
      do i=1,nlats
         if ( latvalues(i) .lt. tmpreal ) then
            if ( i .eq. 1 ) then
               istart(2) = 1
            else
               istart(2) = i
            end if
            exit
         end if
      end do
      tmpreal = minlat
      do i=istart(2),nlats
         if ( latvalues(i) .lt. tmpreal ) then
            icount(2) = i - istart(2)
            exit
         end if
      end do
      if ( icount(2) .eq. 1 ) then
         print *, 'Error in rdhrvegmap, icount(2)=1'
         ierror = -1
         stop 1
      end if
      

      hrmapres = abs(xres1)     
      scaling = nint( xres / hrmapres )
      scaling2 = dble(scaling) * dble(scaling)
      hrmapnlon=icount(1)
      hrmapnlat=icount(2)
      
      !print *,hrmapnlon,hrmapnlat
      !print *,lonvalues(istart(1)),lonvalues(istart(1) + icount(1)-1)

      ! fill global vars hrmaplonvalues, hrmaplatvalues
      allocate(hrmaplonvalues(1:hrmapnlon)) 
      allocate(hrmaplatvalues(1:hrmapnlat)) 
      hrmaplonvalues(:) = 0.
      hrmaplatvalues(:) = 0.
      do i = istart(1), istart(1) + icount(1) - 1
         hrmaplonvalues(i-istart(1)+1) = lonvalues(i)
      end do
      do i = istart(2), istart(2) + icount(2) - 1
         hrmaplatvalues(i-istart(2)+1) = latvalues(i)
      end do

! Read data
! ----------
      ierr = NF_INQ_VARID(idies,varname,idvar)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in rdhrvegmap, cannot get var id'
         print *, NF_STRERROR(ierr)
         ierror = -1
         stop 1
      end if
      allocate(values(1:icount(1),1:icount(2))) 
      ierr = NF_GET_VARA_INT(idies,idvar,istart,icount,values)
      if (ierr .ne. NF_NOERR) then 
         print *, 'Error in rdhrvegmap, cannot get var values'
         print *, NF_STRERROR(ierr)
         ierror = -1
         stop 1
      end if


! Process data
! ----------
     ! allocate vars
      allocate(vegtypecount(npoi1,nvegtype+1)) ! 15 vegtypes + others 
      vegtypecount(:,:) = 0
      allocate(tmpvegtypecount(nvegtype)) 
      tmpvegtypecount(:) = 0
      allocate(xinvegbak(npoi)) 
      xinvegbak = xinveg
      xinveg(:) = 0
      allocate(tilefrac1(npoi))
      allocate(tilefrac2(npoi1))
      tilefrac(:) = 0.
      tilefrac1(:) = 0.
      tilefrac2(:) = 0.
      waterfrac(:) = 0.
      ihrtileparentoffset = 32768
      allocate(ihrtileparent(hrmapnlat,hrmapnlon)) ! global
      ihrtileparent(:,:) = - ihrtileparentoffset
      !allocate(ihrtileparent(hrmapnlon*hrmapnlat)) ! global
      !ihrtileparent(:) = -ihrtileparentoffset

      ! count vegtype total for each model grid point
      do inpoi = 1,npoi1

         ! TODO make this equation simpler... it also seems fragile
         ! could probably just loop over hrmaponvalues/hrmaplatvalues
         jb = nint(( latvalues(istart(2)) - latscale(latindex(inpoi)) - (yres/2.0) ) / yres ) * scaling
         do j=jb+1,jb+scaling

            ib = nint(( lonscale(lonindex(inpoi)) - lonvalues(istart(1))  - (xres/2.0) ) / xres ) * scaling
            do i=ib+1,ib+scaling
               tmpint = values(i,j)
               ! transform some values to INLAND values
               ! wetlands = water
               if (tmpint .eq. 16) then
                  tmpint = watervegtype
                  values(i,j) = watervegtype
               end if
               ! water goes to waterfrac
               if (tmpint .eq. watervegtype) then
                  waterfrac(inpoi) = waterfrac(inpoi) + 1
               ! > nvegtype goes to nvegtype+1
               else if (tmpint .gt. nvegtype .or. tmpint .lt. 1) then
                  !print *,i,j,tmpint,'out of nvegtype range'
                  tmpcount = tmpcount+1
                  vegtypecount(inpoi,nvegtype+1) = vegtypecount(inpoi,nvegtype+1) + 1
               else
                  vegtypecount(inpoi,tmpint) = vegtypecount(inpoi,tmpint) + 1
                  !ihrtileparent(i,j) = inpoi - ihrtileparentoffset
               end if

            end do ! i

         end do ! j

      end do ! inpoi

      if ( env_debug .gt. 3 ) then
         do inpoi=1,npoi1
            print *,"1",inpoi,sum(vegtypecount(inpoi,:)),dble(sum(vegtypecount(inpoi,:)))/scaling2,waterfrac(inpoi),waterfrac(inpoi)/scaling2
            print *,inpoi,vegtypecount(inpoi,:)
            !print *,inpoi,maxloc(vegtypecount(inpoi,1:nvegtype)),maxval(vegtypecount(inpoi,1:nvegtype))
         end do
      end if

      ! find dominant and other tile categories, write to xinveg
      do inpoi=1,npoi1
         tmpvegtypecount = vegtypecount(inpoi,1:nvegtype)
         domtile = 0

         do ilpt = 1,mlpt
            !print *,ilpt,sum(tmpvegtypecount)
            if ( sum(tmpvegtypecount) .ne. 0 ) then
               max_loc = maxloc(tmpvegtypecount) ! get current max value
               if ( ilpt .eq. 1 ) domtile = max_loc(1)              
               isubgrid = subgrid_get_index(inpoi,ilpt)
               xinveg( isubgrid ) = max_loc(1)
               xinveg( isubgrid ) = max_loc(1)
               tilefrac1( isubgrid ) = real(tmpvegtypecount(max_loc(1))) / scaling2
               tmpvegtypecount( max_loc(1) ) = 0 ! zero this element to get next max value
            ! if sum == 0 then we got all vegtypes for this grid box
            else
               ! if this occurs with first subgrid tile, it means we got no valid vegtypes,
               ! fill with global vegtype map
               if ( ilpt .eq. 1 ) then 
                  domtile = xinvegbak( subgrid_get_index(inpoi,1) )
                  if ( env_debug .gt. 1 ) then
                     print *, 'FIXME: rdhrvegmap, cannot get dominant tile for npoi=',inpoi,&
                          lonscale(lonindex(inpoi)),latscale(latindex(inpoi))
                     print *, '       falling back to global map, using vegtype=',domtile
                  end if
                  isubgrid = subgrid_get_index(inpoi,1)
                  xinveg( isubgrid ) = domtile
                  tilefrac( isubgrid ) = 1.
                  tilefrac1( isubgrid ) = 1.
                  do i = 2,mlpt
                     isubgrid = subgrid_get_index(inpoi,i)
                     xinveg( isubgrid ) = domtile
                     tilefrac( isubgrid ) = 0.
                     !tilefrac1( isubgrid ) = 0.
                  end do
               ! else, fill with domtile value
               else
                  xinveg( subgrid_get_index(inpoi,ilpt) ) = domtile
               end if
            end if
         end do !ilpt
         !update tilefrac2 element with out-of bounds ids+water
         tilefrac2( inpoi ) = vegtypecount(inpoi,nvegtype+1) + waterfrac( inpoi )
      end do !inpoi

     ! fill out-of-range vegtype values with dominant id
     ! fill values that are not one of the tiles with dominant id
      
      do inpoi = 1,npoi1
         domtile = xinveg(inpoi)
         !print *,'inpoi=',inpoi,domtile,lonscale(lonindex(inpoi)),latscale(latindex(inpoi))
         jb = nint(( latvalues(istart(2)) - latscale(latindex(inpoi)) - (yres/2.0) ) / yres ) * scaling
         do j=jb+1,jb+scaling
            ib = nint(( lonscale(lonindex(inpoi)) - lonvalues(istart(1))  - (xres/2.0) ) / xres ) * scaling
            do i=ib+1,ib+scaling
               ! set hr tile parent as dominant tile by default
               ihrtileparent(j,i) = inpoi - ihrtileparentoffset
               tmpint = values(i,j)
               ! illegal values
               ! TODO transform some values to INLAND values
               !if (tmpint .gt. nvegtype .or. tmpint .lt. 1) then
               !   if ( tmpint .eq. 19 ) then
               !      tmpint = 
               !   if 
               !end if
               ! water pixel, do nothing
               if (tmpint .eq. watervegtype) then
                  ihrtileparent(j,i) = - ihrtileparentoffset ! or dominant?
               ! out of range pixel
               else if ( tmpint .gt. nvegtype .or. tmpint .lt. 1) then
                  values(i,j) = domtile
                  vegtypecount(inpoi,nvegtype+1) = vegtypecount(inpoi,nvegtype+1) - 1
                  vegtypecount(inpoi,domtile) = vegtypecount(inpoi,domtile) + 1
               else
                  ! test if value is one of the tiles, if not add to dominant tiles
                  ok = .FALSE.
                  do ilpt = 1,mlpt
                     isubgrid = subgrid_get_index(inpoi,ilpt)
                     if (tmpint .eq. xinveg( isubgrid )) then
                        ok = .TRUE.
                        ! also set hr tile parent as this subtrid tile
                        ihrtileparent(j,i) = isubgrid - ihrtileparentoffset
                        ! leave do loop as we found it
                        exit
                     end if
                  end do
                  if ( .not. ok ) then
                     !print *,i,j,tmpint,'non-tile'
                     vegtypecount(inpoi,tmpint) = vegtypecount(inpoi,tmpint) - 1
                     vegtypecount(inpoi,domtile) = vegtypecount(inpoi,domtile) + 1
                     tilefrac2( inpoi ) = tilefrac2( inpoi ) + 1 
                  end if
               end if
            end do
         end do
      end do

      if ( env_debug .gt. 3 ) then
         do inpoi=1,npoi1
            !print *,inpoi,vegtypecount(inpoi,:),waterfrac(inpoi),sum(vegtypecount(inpoi,:))
            print *,"2",inpoi,sum(vegtypecount(inpoi,:)),dble(sum(vegtypecount(inpoi,:)))/scaling2,waterfrac(inpoi),waterfrac(inpoi)/scaling2
            print *,inpoi,vegtypecount(inpoi,:)
            !print *,inpoi,maxloc(vegtypecount(inpoi,1:nvegtype)),maxval(vegtypecount(inpoi,1:nvegtype))
         end do
      end if

      ! update waterfrac and tilefrac
      waterfrac = waterfrac / scaling2
      tilefrac2 = tilefrac2 / scaling2
      do inpoi=1,npoi1
         tmpvegtypecount = vegtypecount(inpoi,1:nvegtype)
         do ilpt = 1,mlpt
            isubgrid = subgrid_get_index(inpoi,ilpt)
            if ( sum(tmpvegtypecount) .ne. 0 ) then
               max_loc = maxloc(tmpvegtypecount) ! get current max value
               tilefrac( isubgrid ) = real(tmpvegtypecount(max_loc(1))) / scaling2
               tmpvegtypecount( max_loc(1) ) = 0 ! zero this element to get next max value
            else
               tilefrac( isubgrid ) = 0.
            end if
            if ( env_debug .gt. 3 ) then
              print *,inpoi,ilpt,max_loc(1),tilefrac(isubgrid),tilefrac1(isubgrid),tilefrac2(inpoi)
            end if
         end do
      end do


! Close file
! ----------
     ierr = NF_CLOSE(idies)
     if (ierr .ne. NF_NOERR) then
        print *, 'Error in rdhrvegmap'
        print *, NF_STRERROR(ierr)
        ierror = -1
        stop 1
     end if


! De-allocate local memory
! ----------
     if ( allocated(lonvalues) ) deallocate(lonvalues)
     if ( allocated(latvalues) ) deallocate(latvalues)
     if ( allocated(values) ) deallocate(values)
     if ( allocated(vegtypecount) ) deallocate(vegtypecount)
     if ( allocated(tmpvegtypecount) ) deallocate(tmpvegtypecount)
     if ( allocated(xinvegbak) ) deallocate(xinvegbak)

     ! TMP TMP
     !if ( allocated(ihrtileparent) ) deallocate(ihrtileparent)

     print *,'INFO: done processing high-res map'
     return

end subroutine rdhrvegmap

#endif
