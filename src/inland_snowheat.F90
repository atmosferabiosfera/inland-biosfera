#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine snowheat (tlsub, fhtop, sflo, xl, chl, kpti, kptj)
! ---------------------------------------------------------------------
! sets up call to tridiaib to solve implicit snow heat conduction,
! using snow temperatures in tsno (in comsno). adds an extra
! buried-lower-veg layer to the bottom of the snow with 
! conduction coefficient conbur/xl and heat capacity chl*xl
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comsno
      use inland_comsoi

      implicit none
!-----------------------------------------------------------------------
! input-output variables
      integer kpti         ! index of 1st point of little vector 
                           ! in big lpt vector
      integer kptj         ! index of last point of little vector
 
      real*8 chl           ! specific heat of lower veg per l/s area (supplied)

      real*8 tlsub(lbeg:lend), & ! temperature of buried lower veg (supplied, returned)
             fhtop(kpti:kptj), & ! heat flux into top snow layer from atmos (supplied)
             sflo(kpti:kptj,nsnolay+2), & ! downward heat flow across layer boundaries (returned)
             xl(kpti:kptj)                ! (lai(i,1)+sai(i,1))*fl(i), lower-veg density(supplied)

! local variables
      integer npt        ! number of points in little vector
      integer k, i, j, & ! loop indices
              km1,     & ! used to avoid layer 0
              kp1        ! used to avoid layer nsnolay+2
!             
      real*8 rimp,   &   ! implicit fraction of the calculation (0 to 1)
             conbur, &   ! conduction coeff of buried lower veg layer for unit density xl=(lai+sai)*fl, in w m-2 k-1
             hfake,  &   ! arbitrary small thickness to allow processing for zero snow. (doesn't use index since tridia
                         ! not set up for index.)
             rwork,  &   ! to compute matrix diagonals and right-hand side
             dt,     &   ! '
             dti         ! '

      real*8 con(mpt,nsnolay+2),  & ! conduction coefficents between layers
             temp(mpt,nsnolay+1), & ! combined snow and buried-veg temperatures
             d1(mpt,nsnolay+1),   & ! diagonals of tridiagonal systems of equations 
             d2(mpt,nsnolay+1),   & !  diagonals of tridiagonal systems of equations 
             d3(mpt,nsnolay+1),   & ! diagonals of tridiagonal systems of equations 
             rhs(mpt,nsnolay+1),  & ! right-hand sides of systems of equations
             w1(mpt,nsnolay+1),   & ! work array needed by tridia
             w2(mpt,nsnolay+1)

! conbur (for xl=1) is chosen to be equiv to 10 cm of snow
      data rimp, conbur, hfake /1.0, 2.0, .01/
      npt = kptj - kpti + 1

! initialize variables used by tridiaib
      d1(:,:)=0.
      d2(:,:)=0.
      d3(:,:)=0.
      rhs(:,:)=0.
      w1(:,:)=0.
      w2(:,:)=0.

! copy snow and buried-lower-veg temperatures into combined
! array temp
      call scopya(npt, tlsub(kpti), temp(1,nsnolay+1))

      do k = 1, nsnolay
         do i = kpti, kptj
            j = i - kpti +1
            temp(j,k) = tsno(i,k)
         end do
      end do

      con(:,:)=0.
! set conduction coefficients between layers
! it is -between- layers, i.e. the first layer is excluded as we calculate 
! the conductancy between the current and previous layers.
      do 100 k=2,nsnolay+2
         if (k.le.nsnolay) then
            rwork = 0.5 / consno
            do 102 i = kpti, kptj
               j = i - kpti +1
               con(j,k) = 1. / (max(hsno(i,k-1),hfake)*rwork + max(hsno(i,k), &
                          hfake)*rwork )
102         continue
         elseif (k.eq.nsnolay+1) then
            rwork = 0.5 / consno
            do 104 i = kpti, kptj
               j = i - kpti +1
               con(j,k) = 1. / (max(hsno(i,k-1),hfake)*rwork + &
                          0.5*xl(i)/conbur )
104         continue
         elseif (k.eq.nsnolay+2) then
            rwork = 0.5 / conbur
            do 106 i = kpti, kptj
               j = i - kpti +1
               con(j,k) = 1. / (   xl(i)*rwork + 0.5*hsoi(1) / consoi(i,1) )
106         continue
         endif
100   continue

! set matrix diagonals and right-hand side. for layer nsnolay+1
! (buried-lower-veg layer), use explicit contact with soil, and
! multiply eqn through by xl*chl/dtime to allow zero xl.
      do 200 k=1,nsnolay+1
         km1 = max (k-1,1)
         kp1 = min (k+1,nsnolay+1)
         if (k.le.nsnolay) then
            rwork = dtime /(rhos*cice)
            do 202 i = kpti, kptj
               j = i - kpti +1
               dt = rwork / (max(hsno(i,k),hfake))
               d1(j,k) =    - dt*rimp* con(j,k)
               d2(j,k) = 1. + dt*rimp*(con(j,k)+con(j,k+1))
               d3(j,k) =    - dt*rimp* con(j,k+1)
               rhs(j,k)=temp(j,k)+dt*((1.-rimp)*con(j,k)*(temp(j,km1)- &
                        temp(j,k))+(1.-rimp)*con(j,k+1)*(temp(j,kp1)-  &
                        temp(j,k)))
202         continue
            if (k.eq.1) then 
               rwork = dtime /(rhos*cice)
               do 204 i = kpti, kptj
                  j = i - kpti +1
                  dt = rwork / (max(hsno(i,k),hfake))
                  rhs(j,k) = rhs(j,k) + dt*fhtop(i)
204            continue
            endif
         elseif (k.eq.nsnolay+1) then
            rwork = chl / dtime
            do 206 i = kpti, kptj
               j = i - kpti +1
               dti = xl(i)*rwork
               d1(j,k) =     -  rimp* con(j,k)
               d2(j,k) = dti +  rimp*(con(j,k)+con(j,k+1))
               d3(j,k) = 0.
               rhs(j,k) = dti*temp(j,k)+((1.-rimp)*con(j,k)*(temp(j,km1)- &
                          temp(j,k)) + con(j,k+1)*(tsoi(i,1)-             &
                          (1.-rimp)*temp(j,k)) )
206         continue
         endif
200   continue

! solve the tridiagonal systems
      call tridiaib (npt, mpt,nsnolay+1,d1,d2,d3,rhs, temp, w1,w2)

! deduce downward heat fluxes between layers
      do i = kpti, kptj
         sflo(i,1) = fhtop(i)
      end do

      do 400 k=1,nsnolay+1
         if (k.le.nsnolay) then
            rwork = rhos*cice/dtime
            do 402 i = kpti, kptj
               j = i - kpti +1
               sflo(i,k+1) = sflo(i,k)-rwork*hsno(i,k)*(temp(j,k)-tsno(i,k))
402         continue
         else
            rwork = chl/dtime
            do 404 i = kpti, kptj
               j = i - kpti +1
               sflo(i,k+1)=sflo(i,k)-xl(i)*rwork*(temp(j,nsnolay+1)-tlsub(i))
404            continue
         endif
400   continue

! copy temperature solution to tsno and tlsub, but not for
! points with no snow
      do 500 k=1,nsnolay
         do 502 i =  kpti, kptj
            j = i - kpti +1
            if (fi(i).gt.0.) tsno(i,k) = temp(j,k) 
502      continue
500   continue

      do 510 i = kpti, kptj
         j = i - kpti +1
         if (fi(i).gt.0.) tlsub(i) = temp(j,nsnolay+1)
510   continue
      return
end subroutine snowheat
