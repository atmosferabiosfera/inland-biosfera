#include "inland_config.h"
! ---------------------------------------------------------------------
subroutine soilheat (otsoi, owsoi, c0pud, fhtop, c1pud, &
                     kpti, kptj)
! ---------------------------------------------------------------------
!
!        sets up call to tridiaib to solve implicit soil/ice heat 
!        conduction, using layer temperatures in tsoi (in comsoi).
!        the heat flux due to liquid flow previously calculated
!        in soilh2o is accounted for. lower bc is conductive flux = 0
!        for soil (although the flux due to liquid drainage flow can
!        be > 0)
!
! ---------------------------------------------------------------------
      use inland_parameters
      use inland_comsoi

      implicit none
!-----------------------------------------------------------------------
! input-output variables
      integer kpti  ! index of 1st point of little vector in big lpt vector
      integer kptj  ! index of last point of little vector

! otsoi: soil/ice temperatures at start of timestep (redundant with tsoi, but 
!       passed to be consistent with soilh2o)  
! owsoi: soil moistures at start of timestep (before soilh2o)
! c0pud: layer heat capacity due to puddles (=0 except for top)
! c1pud: updated c0pud
! fhtop: heat flux into top layer from atmos
! wflo: downward h2o  flow across layer boundaries
      real*8 otsoi(lbeg:lend,nsoilay), owsoi(lbeg:lend,nsoilay), &   
             c0pud(lbeg:lend,nsoilay), c1pud(lbeg:lend,nsoilay), &
             fhtop(lbeg:lend)  

! local variables
      integer npt               ! number of points in little vector
      integer k, i, km1, kp1, j

      real*8 rimp,           &  ! implicit fraction of the calculation (0 to 1)
             rwork, rwork1,  &  ! work variables
             rwork2, rwork3, &  ! work variables
             t

! whflo: downward heat fluxes across layer bdries due to h2o movement calculated in soilh2o
! con: conduction coeffs between layers
! c0: specific heats at start of timestep
! c1: specific heats at end of timestep
! d1,d2,d3: diagonals of tridiagonal systems of equations
! rhs: right-hand sides of systems of equations
! tsoiw: translation of tsoi(kpti:kptj) to tsoiw(1:npt)
! w1,w2: work arrays needed by tridia
      real*8 whflo(lbeg:lend,nsoilay+1), con(lbeg:lend,nsoilay+1), &
       c0(lbeg:lend,nsoilay), c1(lbeg:lend,nsoilay), d1(mpt,nsoilay), &
       d2(mpt,nsoilay), d3(mpt,nsoilay), rhs(mpt,nsoilay), tsoiw(mpt,nsoilay), &
       w1(mpt,nsoilay), w2(mpt,nsoilay)

      data rimp /1.0/
!-----------------------------------------------------------------------
! set conduction coefficient between layers, and heat fluxes
! due to liquid transport
!
! top layer
      k = 1
      do 100 i = kpti, kptj
         con(i,k) = 0.0
         whflo(i,k) = wflo(i,k) * ch2o * tsoi(i,k)
100   continue

! middle layers
      do 110 k = 2, nsoilay
         do 120 i = kpti, kptj
            con(i,k) = 1. / (0.5*(hsoi(k-1)/consoi(i,k-1)+hsoi(k)/consoi(i,k)))
            t=(hsoi(k)*tsoi(i,k-1) + hsoi(k-1)*tsoi(i,k))/(hsoi(k-1) + hsoi(k))
            whflo(i,k) = wflo(i,k) * ch2o * t
120      continue
110   continue

! bottom layer
      k = nsoilay + 1
      do 130 i = kpti, kptj
         con(i,k) = 0.0
         whflo(i,k) = wflo(i,k) * ch2o * tsoi(i,k-1)
130   continue

! set diagonals of matrix and right-hand side. use old and
! new heat capacities c0, c1 consistently with moisture fluxes
! whflo computed above, to conserve heat associated with 
! changing h2o amounts in each layer
      do 200 k = 1, nsoilay
         km1 = max (k-1,1)
         kp1 = min (k+1,nsoilay)
         do 210 i= kpti, kptj
            j = i - kpti + 1
            rwork1 = (1.-poros(i,k))*csoi(i,k)*rhosoi(i,k)
            rwork2 = poros(i,k)*(1.-wisoi(i,k))*ch2o*rhow
            rwork3 = poros(i,k)*wisoi(i,k)*cice*rhow
            c0(i,k) = c0pud(i,k) + (rwork1 + rwork2*owsoi(i,k) + rwork3)*hsoi(k)
            c1(i,k) = c1pud(i,k) + (rwork1 + rwork2*wsoi(i,k) + rwork3)*hsoi(k)
            rwork = dtime/c1(i,k)
            d1(j,k) =    - rwork * rimp * con(i,k)
            d2(j,k) = 1. + rwork * rimp * (con(i,k)+con(i,k+1))
            d3(j,k) =    - rwork * rimp * con(i,k+1)
            rhs(j,k) = (c0(i,k)/c1(i,k))*tsoi(i,k) + rwork &
                       * ( (1.-rimp)*con(i,k)  *(tsoi(i,km1)-tsoi(i,k)) &
                       + (1.-rimp)*con(i,k+1)*(tsoi(i,kp1)-tsoi(i,k)) &
                       + whflo(i,k) - whflo(i,k+1) )
210      continue
         if (k.eq.1) then
            do 220 i= kpti, kptj
               j = i - kpti + 1
               rhs(j,k) = rhs(j,k) + (dtime/c1(i,k))*fhtop(i)
220         continue
         endif
200   continue

! solve systems of equations
      npt = kptj - kpti + 1
      call tridiaib(npt,npt,nsoilay,d1,d2,d3,rhs,tsoiw,w1,w2)

! Translates tsoiw(1 --> npt) from tridiaib to kpti --> kptj 
      do k = 1, nsoilay
         do i = kpti, kptj
            j = i - kpti + 1
            tsoi(i,k) = tsoiw(j,k)
         end do
      end do

! deduce downward heat fluxes between layers
      do i = kpti, kptj
         hflo(i,1) = fhtop(i)
      end do

      do 300 k = 1, nsoilay
         do 310 i = kpti, kptj
            hflo(i,k+1)=hflo(i,k)-(c1(i,k)*tsoi(i,k)-c0(i,k)*otsoi(i,k))/dtime
310      continue
300   continue
      return
end subroutine soilheat
