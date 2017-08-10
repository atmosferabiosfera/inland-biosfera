#include "inland_config.h"
! ---------------------------------------------------------------------
      subroutine diaginit(idiag, diaglat, diaglon)
! ---------------------------------------------------------------------
! CALLED BY inidiag
! CALLS nothing
! Initialize diagnostic output files
!
      use inland_compar.F90
      use inland_comwork.F90
      use inland_comdiag.F90
!
! local variables
!
      character*10 filen          ! name of diagnostic file
      character*10 varname(nvars) ! variable names
!
      character*10 header(16,nfiles)  ! column headers for output file
!
      integer idiag,  &            ! number of diganostic files
              ifile,   &           ! counter
              ivar,     &          ! counter
              kolumn,    &         ! column in output file
              lat(nfiles),&
              lon(nfiles),&
              countvar
!
      real    alat0,&
              alon0
!
      real    diaglat(nfiles),&      ! latitude of diagnostic point
              diaglon(nfiles)        ! longitude of diagnostic point
!
      data varname /'adco2mic','adco2root','adco2soi','adneetot',
                    'adtsoic','adwsoic','amco2root','amco2soi',
                    'ancl3','ancl4','ancub','ancuc','asurd1',
                    'asurd2','asuri1','asuri2','ayco2soi',
                    'biomass1','biomass2','biomass3','biomass4',
                    'biomass5','biomass6','biomass7','biomass8',
                    'biomass9','biomass10','biomass11','biomass12',
                    'cloud','coszen','fi','fira','firb','frac1','frac2',
                    'frac3','frac4','frac5','frac6','frac7','frac8',
                    'frac9','frac10','frac11','frac12','gadjust',
                    'gdrain','ginvap','grunof','gsuvap','gtrans',
                    'gtransl','gtransu','hsno1','hsno2','hsno3',
                    'hsnotop','plai1','plai2','plai3','plai4','plai5',
                    'plai6','plai7','plai8','plai9','plai10','plai11',
                    'plai12','precip','psurf','qa','raina','snowa',
                    'solad1','solad2','solai1','solai2','ta','td','tl',
                    'totalit','totcmic','totcsoi','totrlit','trng',
                    'ts','tsno1','tsno2','tsno3','tsoi1','tsoi2',
                    'tsoi3','tsoi4','tsoi5','tsoi6','tu',
                    'wet','wipud','wisoi1','wisoi2','wisoi3','wisoi4',
                    'wisoi5','wisoi6','wliql',
                    'wliqs','wliqu',
                    'wpud','wsnol','wsnos','wsnou','wsoi1','wsoi2',
                    'wsoi3','wsoi4','wsoi5','wsoi6',
                    'biomass13', 'biomass14', 'biomass15',
                    'frac13', 'frac14', 'frac15',
                    'plai13', 'plai14', 'plai15' /
!
! Initialize diagnostic output files
!
      do 10 ifile = 1,nfiles
        kolumn = 4
        countvar = 0
        do 20 ivar = 1,nvars
          if ( ldiag(ivar,ifile) .eq. 1 ) then
            kolumn = kolumn + 1
            header(kolumn,ifile) = varname(ivar)
            countvar = countvar + 1
          end if
 20     continue
          if (countvar .lt. 12) then
            do 25 i = countvar+5,16
              header(i,ifile) = 'empty'
 25         continue
          end if
 10   continue
!    
      do 30 i = 1, idiag
        filen = 'ibis.diag'
        write(filen(10:10),'(i1)')i-1
        open (i+20, status='unknown', file=filen, access='append')
        write (i+20,3000) diaglat(i), diaglon(i)
        write (i+20,3100) diagstart(i), diagend(i)
        write (i+20,3200) (header(kolumn,i),kolumn=5,16)
!
 3000 format ('%Cell_latitude ', f8.2,1x,' Cell_longitude ',f8.2)
 3100 format ('%Begin_year ',i8,' End_year ',i8)
 3200 format ('%Year',2x,'Month',4x,'Day',4x,'Step',3x,12(a12,1x))
!
 30   continue
!
! Determine the grid cells corresponding to the chosen coordinates
!
      alat0 = latscale(1)
      alon0 = lonscale(1)
!
      do 40 i = 1, idiag
        lat(i) = 1 + int( (alat0 - diaglat(i))/yres + 0.5 )
        lat(i) = min(lat(i),nlatsub)
        lat(i) = max(lat(i),1)
        lon(i) = 1 + int( (diaglon(i) - alon0)/xres + 0.5 )
        lon(i) = min(lon(i),nlonsub)
        lon(i) = max(lon(i),1)
        ndiagpt(i) = -999.
      do 50 i = lbeg, lend
          if (latindex(n) .eq. lat(i) .and. 
     >        lonindex(n) .eq. lon(i)) then
             ndiagpt(i) = n
          end if
 50    continue
       if (ndiagpt(i) .lt. 1) then
          print *, 'Warning: diag point ',i,':',diaglon(i),diaglat(i),
     >     ' is not a land point.  Ignoring'
       end if
 40   continue
!
      return
      end subroutine diaginit
