! TODO: make inland_comsat.F90 module in place of this ugly patch.
!
! TODO fix all warnings such as
! ../include/inland_comsat.h:10.22:
!    Included at inland_diurnal.F90:265:
!
! real*8 t, p1, e1, tair
!                       1
! Warning: Unused variable 'tair' declared at (1)
!
! ------
! comsat
! ------
!
! ---------------------------------------------------------------------
! statement functions associated with module comsatp parameters
! ---------------------------------------------------------------------
!
real*8 t, p1, e1, tair
#ifndef SETSOI_COMSAT
# if defined(STEPH2O2_COMSAT)
#  define SETSOI_COMSAT
# endif /* STEPH2O2_COMSAT */
#endif /* !SETSOI_COMSAT */

#ifndef WGEN_COMSAT
# if defined(STOMATAIB_COMSAT)
#  define WGEN_COMSAT
# endif /* STOMATAIB_COMSAT */
#endif /* WGEN_COMSAT */

#ifdef WGEN_COMSAT
real*8 tsatl, tsati, esat, qsat
#else /* WGEN_COMSAT */
# ifdef SETSOI_COMSAT
real*8 hvapf, hsubf
# else /* SETSOI_COMSAT */
#  ifdef TWET3_COMSAT
real*8 tsatl, tsati, esat, desat, qsat, hvapf
#  else /* TWET3_COMSAT */
real*8 tsatl, tsati, esat, desat, qsat, dqsat, q1, hvapf, hsubf
#  endif /* TWET3_COMSAT */
# endif /* SETSOI_COMSAT */
#endif /* WGEN_COMSAT */

#ifndef SETSOI_COMSAT
! statement functions tsatl,tsati are used below so that lowe's
! polyomial for liquid is used if t gt 273.16, or for ice if 
! t lt 273.16. also impose range of validity for lowe's polys.
tsatl(t) = min (dble(100.), max (t-273.16, dble(0.)))
tsati(t) = max (dble(-60.), min (t-273.16, dble(0.)))

! statement function esat is svp in n/m**2, with t in deg k. 
! (100 * lowe's poly since 1 mb = 100 n/m**2.)
esat(t) = dble(100.)*(merge(asat0, bsat0, t.ge.273.16) +          &
                 tsatl(t)*(asat1 + tsatl(t)*(asat2 + tsatl(t) *   &
          (asat3 + tsatl(t)*(asat4 + tsatl(t)*(asat5 + tsatl(t) * &
          asat6))))) +                                            &
          tsati(t)*(bsat1 + tsati(t)*(bsat2 + tsati(t) *          &
          (bsat3 + tsati(t)*(bsat4 + tsati(t)*(bsat5 + tsati(t) * &
          bsat6))))))

! statement function qsat is saturation specific humidity,
! with svp e and ambient pressure p in n/m**2. impose an upper
! limit of 1 to avoid spurious values for very high svp
! and/or small p
qsat(e1, p1) = 0.622*e1/max(p1-(1.0-0.622)*e1, 0.622*e1)
#endif /* SETSOI_COMSAT */

#if !defined WGEN_COMSAT
# if !defined SETSOI_COMSAT && !defined TWET3_COMSAT
! statement function dsat is d(sat. spec. humidity)/dt, with t 
! in deg k, and neglecting q in denominator of the q(vp) reln.
! (100 * lowe's poly since 1 mb = 100 n/m**2.)
desat(t) = 100.*(merge(csat0, dsat0, t.ge.273.16) +                       &
                 tsatl(t)*(csat1 + tsatl(t)*(csat2 + tsatl(t) *          &
           (csat3 + tsatl(t)*(csat4 + tsatl(t)*(csat5 + tsatl(t) * &
                                               csat6))))) +               &
                  tsati(t)*(dsat1 + tsati(t)*(dsat2 + tsati(t) *          &
           (dsat3 + tsati(t)*(dsat4 + tsati(t)*(dsat5 + tsati(t)* &
                                                dsat6))))))

! statement function dqsat is d(qsat)/dt, with t in deg k and q
! in kg/kg (q is *saturation* specific humidity)
dqsat(t, q1) = desat(t)*q1*(1. + q1*(1./0.622-1.)) / esat(t)
# endif /* !SETSOI_COMSAT && !TWET3_COMSAT */
!
! statement functions hvapf, hsubf correct the latent heats of
! vaporization (liquid-vapor) and sublimation (ice-vapor) to
! allow for the concept that the phase change takes place at
! 273.16, and the various phases are cooled/heated to that 
! temperature before/after the change. this concept is not
! physical but is needed to balance the "black-box" energy 
! budget. similar correction is applied in convad in the agcm
! for precip. needs common comgrd for the physical constants.
! argument t is the temp of the liquid or ice, and tair is the
! temp of the delivered or received vapor.
!
hvapf(t,tair) = hvap + cvap*(tair-273.16) - ch2o*(t-273.16)

# ifndef TWET3_COMSAT
hsubf(t,tair) = hsub + cvap*(tair-273.16) - cice*(t-273.16)
!
# endif /* TWET3_COMSAT */
#endif /* WGEN_COMSAT */
