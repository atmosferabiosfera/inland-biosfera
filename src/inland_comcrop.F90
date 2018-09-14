#include "inland_config.h"
module inland_comcrop

      implicit none
      public
      save
!
!--------
! comcrop last edited by C. Kucharik 01.30.03
!--------
!
! TODO replace isimagro tests with point-specific tests on landusetype (except in main)
!      this is necessary to mix natural and agro vegetation in different subtrids
!      - Pousa's job!!!
! also - modify params/ reading and soil levels to be able to mix agro and natural
      integer   isimagro,        &  ! 0: no agro crops, 1 agro crops - define either icroptype or cropsfile
                icroptype,       &  ! define crop type to use in all grid points - 0: no crops / 13: soybean / 14: maize / 15: wheat / 16: sgc / 17: oil palm
! TODO iwheattype should be a grid-point variable instead...
                iwheattype,      &  ! 1: spring wheat 2: winter wheat - used if icroptype = wheat
                irotation,       &  ! 0: no crops rotated (monoculture) 2: two crops rotated 3: 3 crops rotated (see crops.f)
                isoilay,         &  ! soil layer for which drainage and leaching data is output
                nratoon,         &  ! number of ratoon crops
                qgalho,          &  !
                nplant,          &  !
                istyear,         &  ! year to start using daily station data for climate
                istend              ! ending year of reading in daily station data
!
      real*8, dimension(:), allocatable :: consdays, maxcons, gsdays, xirrig, xirriga,&
                                           totirrig, df, vf, cumvd, hdidx, gddfzcorn, &
                                           gddfzsgc, gddfzsoy, gddfzwht, conspdate,   &
                                           conshybrid, ncyears, cropy, ik, gddpl15,   &
                                           gddemerg, rm, croptype

!
!            consdays(npoi),        ! counter keeping track of number of days between freezing events
!            maxcons(npoi),         ! maximum number of consecutive days between freezing events
!            gsdays(npoi),          ! number of days in the growing season
!            xirrig(npoi),          ! irrigated water application rate (mm/day) to crops
!            xirriga(npoi),         ! irrigated application rate per timestep
!            totirrig(npoi),        ! annual total irrigation applied (mm/yr)
!            df(npoi),              ! photoperiod factor for wheat crops
!            vf(npoi),              ! vernalization factor for wheat crops
!            cumvd(npoi),           ! cumulative vernalization days
!            hdidx(npoi),           ! cold-hardening index for wheat (varies from 0-2: stage 1 : 0-1, stage 2: 1-2)
!            gddfzcorn(npoi),       ! seasonal gdd (base 8 C) in between last freeze and first freeze
!            gddfzsgc(npoi),        ! seasonal gdd (base 12 C) in between last freeze and first freeze
!            gddfzsoy(npoi),        ! seasonal gdd (base 10 C) in between last freeze and first freeze
!            gddfzwht(npoi),        ! seasonal gdd (base 0 C) in between last freeze and first freeze
!            conspdate(npoi),       ! constant planting date based on specified time period
!            conshybrid(npoi),      ! constant hybrid (gdd) based on specified time period
!            ncyears(npoi),         ! number of the ith crop year
!            cropy(npoi),           ! sugarcane crop year since planting
!            ik(npoi),              ! counter for sugarcane crop one and half year GDD
!            gddpl15(npoi),         ! accumulated growing degree days for 1 and half year for sugarcane
!            gddemerg(npoi),        ! accumulated growing degree days at leaf emergence
!            rm(npoi)               ! sugarcane relative maturity GDD/GDDmaturity
!            icropsum(npoi)         ! index - number of crop types planted in each grid cell - replaced by croptype
!            croptype(npoi)         ! index - crop type (same value as corresponding crop pft code, or 0 if none)
!
      integer, dimension(:), allocatable :: cdays, iniday, endday
!
!            cdays(npoi),           ! days since planting in weather.f
!            iniday(npoi),          ! last day of freeze in spring
!            endday(npoi),          ! first day of freeze in fall
!

      real*8    rootd,          &   !  root decaiment constant, controls de exponential decaimento of root/aerial alocation
                laidcs,         &   ! lai decline start (based on GDD/GDDmaturity)
                mxtsuc,         &   ! maximum effect of temperature on sucrose allocation (function af5)
                af1,            &   !
                af2,            &   ! linear and logarithmical functions related to the stem carbon allocation (astem) - sugarcane
                af3,            &   !
                af4,            &   ! linear and logarithmical functions related to the sucrose allocation (reproductive:arepr) - sugarcane
                af5,            &   !
                af6,            &   ! logar. functions related to the temperature impact on sucrose alloc. (arepr) - sugarcane
                sf1,            &   ! sf1: slope of the line af1 and ipf1: intercept point of line af1
                ipf1,           &   !
                sf3,            &   ! sf3: slope of the line af3 and ipf3: intercept point of line af3
                ipf3,           &   !
                ecf2,           &   ! ecf2: coneffic of the log. functio af2 and ipf2: intercept point of log. func af2
                ipf2,           &   !
                ecf4,           &   ! ecf4: coneffic of the log. functio af4 and ipf4: intercept point of log. func af4
                ipf4,           &   !
                ecf5,           &   ! ecf5: coneffic of the log. functio af5 and tmxf5: maxm temp in the log. func af5
                tf5,            &   !
                wf5,            &   ! wf5 - % of  arepr/astem alloud to change as function of temperature - af5  - (from 0 to 1.0)
                ecf6,           &   ! ecf6: conefficient of the logar. functio af6 and ipf6: intercept point of log. func af6
                ipf6,           &   !
                ecf7,           &   !
                ldf,            &   !
                tmld,           &   ! ldf: leaf decline factor (% of aleaf tha can be reduced); tmld: temp that aleaf*ldf is converted to stalk
                laidc,          &   ! laidc: conefficient of the laidecl
                firecane,       &   ! control if the above ground residue is burn or not after harvest
                alphac,         &   !
                elevin,         &   !
                availn,         &   !
                cnmax,          &   ! maximum allowable residue c/n ratio
                gnmin,          &   !
                smax,           &   ! maximum value stressn can have
                ztopmxsgc,      &   ! maximum height of sugar cane canopy
                ztopmxsoy,      &   ! maximum height of soybean canopy
                ztopmxwht,      &   ! maximum height of wheat canopy
                ztopmxmze           ! maximum height of maize canopy


!      real*8, dimension(:), allocatable ::  ptemp, pmintemp, pmmin, pdmin, pcm, pcd,  &
!                                            aerial, huileaf, huigrain
      real*8, dimension(:), allocatable ::  ptemp, pmintemp, pcm, pcd,  &
                                            aerial, huileaf, huigrain
      real*8, dimension(:), allocatable ::  pmmin_temp, pdmin_temp !gabriel abrahao: Those also have a npoi dimension now, the temp ones are just for reading
      real*8, dimension(:,:), allocatable ::  pmmin, pdmin !gabriel abrahao: Those also have a npoi dimension now

!
!            ptemp(npft),           ! minimum 10 day average temperature for planting (K)
!            pmintemp(npft),        ! minimum 10 day average min temp for planting (K)
!!            pmmin(npft),           ! earliest month to plant (month) DEPRECATED
!!            pdmin(npft),           ! earliest day in earliest month to plant (day) DEPRECATED
!            pmmin(npoi,npft),           ! earliest month to plant (month)
!            pdmin(npoi,npft),           ! earliest day in earliest month to plant (day)
!            pcm(npft),             ! Planting Calendar month (month)
!            pcd(npft)              ! Planting Calendar day  (day)
!            aerial(npft),          ! allocation of carbon to the aerial components (leaf, stem, and reprod)
!            huileaf(npft),  &      ! heat unit index needed to attain leaf emergence after planting
!	     huigrain(npft), &      ! heat unit index needed to reach vegetative maturity
!
      real*8, dimension(:,:), allocatable :: corndop, sgcdop, soydop, whtdop, plmdop, &
                                             gddcorn, gddsgc, gddsgcp, gddsoy, gddwht,&
                                             daygddc, daygddsgc,            &
                                             daygdds

!                                             fertmaize, fertsgc, fertsoy, fertwheat, xinhybrid, xinpdate &
!
!            corndop(npoi,2100),    ! day that corn is planted
!            sgcdop(npoi,2100),     ! day that sugarce is planted
!            plmdop(npoi,2100),     ! day that oil palm is planted
!            soydop(npoi,2100),     ! day that soy  is planted
!            whtdop(npoi,2100)      ! day that wheat is planted
!            gddcorn(npoi,2100),    ! hybrid planted for corn
!            gddsgc(npoi,2100),     ! hybrid planted for sugarcane ratoon, 1 year.
!            gddsgcp(npoi,2),       ! factor to get the hybrid planted for sugarcane plating, 1 and half year.
!            gddsoy(npoi,2100),     ! hybrid planted for soy
!            gddwht(npoi,2100),     ! hybrid planted for wheat
!            fertmaize(npoi,100),   ! historical annual average N-fertilizer applied to maize   crops across US Mississippi basin (1950-00)
!            fertsgc(npoi,100),     ! historical annual average N-fertilizer applied to sugarcane crops across US Mississippi basin (1950-00)
!            fertsoy(npoi,100),     ! historical annual average N-fertilizer applied to soybean crops across US Mississippi Basin (1950-00)
!            fertwheat(npoi,100),   ! historical annual average N-fertilizer applied to wheat   crops across US Mississippi Basin (1950-00)
!            ndepfact(npoi,100),    ! historical annual average N-deposition factor applied to equations in biogeochem.f          (1950-00)
!            daygddc(npoi, 366),    ! gdd accumulation for each particular day of the year (corn)
!            daygddsgc(npoi, 366),  ! gdd accumulation for each particular day of the year (sugarcane)
!            xinpdate(npoi, 100),   ! planting date input for corn
!            daygdds(npoi, 366),    ! gdd accumulation for each particular day of the year (soy)
!            daygddw(npoi, 366),    ! gdd accumulation for each particular day of the year (wheat)
!            xinhybrid(npoi, 100)   ! hybrid average input for corn
!
      real*8, dimension(:,:), allocatable :: cntops, cnrootvec, croplive, grnfraccrop,   &
                                             gddplant, gddtsoi, gddmaturity, thrlai,  &
                                             peaklai, hui, phuf, tlai, templai,       &
                                             harvidx, fnleaf, fnstem, fnroot, fngrain,&
                                             fnplant, tnplant, grainn, cumlvs, idpp,  &
                                             dpgf, cropyld, dmleaf, dmstem, dmresidue,&
                                             dmyield, dmcrop, cropn, cropfixn, nconcl,&
                                             nconcs, nconcr, nconcg, leafout,         &
                                             cropplant, croplaimx, residuen, dmroot,  &
                                             hdate, idppout, pdate, crmclim, crmact,  &
                                             crmplant, ccdays, grainday, fertinput,   &
                                             avehybrid, pstart, laidecl, harvdate,    &
                                             idop, iavepdate
!
!            cntops(npoi,npft),     ! cn ratio of plant residue
!            cnrootvec(npoi,npft),     ! cn ratio of plant roots (renamed to avoid conflict with inland's cnroot)
!            croplive(npoi,npft),   ! 0 crops have been planted and living : 1 crops not living
!            grnfraccrop(npoi,npft),! green fraction of aboveground vegetation for crops
!            gddplant(npoi,npft),   ! accumulated growing degree days past planting date for crop = j
!            gddtsoi(npoi,npft),    ! accumulated growing degree days past planting date based on top layer soil temperature
!            gddmaturity(npoi,npft),! accumulated growing degrees needed for plant to reach both vegetative and physiological maturity
!            thrlai(npoi,npft),     ! lai threshold for crops when senescence begins
!            peaklai(npoi,npft),    ! 0: lai not at maximum allowable, 1: lai at peak value allowed
!            hui(npoi,npft),        ! heat unit index
!            phuf(npoi,npft),       ! previous day's value of the heat unit factor
!            tlai(npoi,npft),
!            templai(npoi,npft),
!            harvidx(npoi,npft),    ! end of year harvest index for crop
!            fnleaf(npoi,npft),     ! current fraction of nitrogen in leaf dry matter
!            fnstem(npoi,npft),     ! current fraction of nitrogen in stem dry matter
!            fnroot(npoi,npft),     ! current fraction of nitrogen in root dry matter
!            fngrain(npoi,npft),    ! current fraction of nitrogen in grain dry matter
!            fnplant(npoi,npft),    ! current fraction of nitrogen in entire plant
!            tnplant(npoi,npft),    ! total nitrogen in plant dry matter
!            grainn(npoi,npft),     ! total nitrogen offtake by crop in grain
!            cumlvs(npoi,npft),     ! total number of leaves emerged
!            idpp(npoi,npft),       ! number of days past planting
!            dpgf(npoi,npft),       ! number of days past grain fill
!            cropyld(npoi,npft),    ! crop yield in t/ha
!            dmleaf(npoi,npft),     ! leaf dry matter in Mg/ha
!            dmstem(npoi,npft),     ! stem dry matter in Mg/ha
!            dmresidue(npoi,npft),  ! aboveground leaf and stem residue dry matter in Mg/ha
!            dmyield(npoi,npft),    ! yield dry matter in Mg/ha
!            dmcrop(npoi,npft),     ! total crop dry matter production in Mg/ha
!            cropn(npoi, npft),     ! nitrogen removed by crop in kg/ha/yr
!            cropfixn(npoi,npft),   ! nitrogen fixation by crop in kg/ha/yr
!            nconcl(npoi,npft),     ! end of season N concentration in leaf (fraction)
!            nconcs(npoi,npft),     ! end of season N concentration in stem (fraction)
!            nconcr(npoi,npft),     ! end of season N concentration in root (fraction)
!            nconcg(npoi,npft),     ! end of season N concentration in grain (fraction)
!            leafout(npoi,npft),    ! gdd accumuation index for leaf emergence in crops
!            cropplant(npoi,npft),  ! index keeping track of whether that crop has been planted during current year
!            croplaimx(npoi,npft),  ! maximum attained lai by crop during growing season
!            residuen(npoi, npft),  ! nitrogen contained in aboveground crop residue
!            dmroot(npoi,npft),     ! fine root dry matter (Mg/ha)
!            hdate(npoi,npft),      ! harvest date (real value)
!            idppout(npoi,npft),      ! harvest date (real value)
!            pdate(npoi,npft),      ! planting date (real value)
!            crmclim(npoi,npft),    ! crop relative maturity rating (CRM) from Pioneer regression relationships
!            crmact(npoi,npft),     ! best crop relative maturity rating (CRM) for that year based on total GDD accumulated
!            crmplant(npoi,npft),   ! crop relative maturity rating (CRM) based on GDD accumulated since planting of that year
!            ccdays(npoi,npft),     ! number of consecutive days with minimum temperatures below killing threshold
!            grainday(npoi,npft),   ! day of year that plant reaches grain fill stage
!            fertinput(npoi,npft),  ! annual fertilizer input for crop (kg-N per m-2)
!            avehybrid(npoi,npft),  ! average hybrid planted (GDD)
!            pstart(npoi,npft),     ! day to start plant, relative to the planting calendar (see params.crp)
!            laidecl(npoi, npft),   ! decline in leaf area for crop
!            harvdate(npoi,npft),   ! day of year that crop pft was harvested
!            idop(npoi,npft),       ! day of year that crop was planted
!            iavepdate(npoi,npft),  ! average planting date over time
!
      real*8, dimension(:,:,:), allocatable :: cropout
!                               cropout(npoi,npft,60), ! crop output var

!
      real*8, dimension(:), allocatable :: arepr, astemi, aleafi, cfrac,       &
                                           baset, mxtmp, tkill, hybgdd, gddmin,       &
                                           lfemerg, grnfill, laicons, allconsl,       &
                                           allconss, laimx, arooti, arootf, aleaff,   &
                                           astemf, declfact, mxgddgf, mxdgfi, mxmat,  &
                                           fleafi, fleaf, cgrain, convfact, maxhi,    &
                                           fyield, fnlfmx, fngrmx, sratio, rratio,    &
                                           fnopt, bfact, grainmoisture
      real*8, dimension(:,:), allocatable :: astem

!
!            arepr(npft),           ! fraction allocation to reproductive organs (grain/fruit) in crops
!            astem(npft),           ! fraction allocation to stem in crops (non leaf/ non grain)
!            astemi(npft),          ! fraction allocation to stem in crops at initial shift to grain
!            aleafi(npft),          ! fraction allocation to leaf in crops at initial shift to grain
!            cfrac(npft),           ! fraction of dry matter production that is carbon
!            baset(npft),           ! base temperature used to accumulate gdd for crops
!            mxtmp(npft),           ! maximum gdd accumulation allowed for each crop pft (per day 0 C)
!            tkill(npft),           ! temperature (K) at which crops are killed due to freeze
!            hybgdd(npft),          ! maximum gdd for specified hybrids
!            gddmin(npft),          ! minimum number of annual gdd needed to plant/grow a crop
!            lfemerg(npft),         ! fraction of annual gdd (to reach physiological mat.) for leaf emergence
!            grnfill(npft),         ! fraction of annual gdd (to reach physiological mat.) for grain fill initiation
!            laicons(npft),         ! lai decline factor constant for crops
!            allconsl(npft),        ! leaf allocation decline scaling factor after grain fill initiation
!            allconss(npft),        ! stem allocation decline scaling factor after grain fill initiation
!            laimx(npft),           ! maximum LAI of each crop allowed
!            arooti(npft),          ! initial allocation of carbon to roots
!            arootf(npft),          ! allocation of carbon to roots at end of growing season
!            aleaff(npft),          ! allocation of carbon to leaves at end of growth cycle
!            astemf(npft),          ! allocation of carbon to stems at end of growth cycle
!            declfact(npft),        ! rate of LAI decline after grain fill inititation (dimensionless factor)
!            mxgddgf(npft),         ! maximum gdd allowed past grain fill inititiation for phys. maturity
!            mxdgfi(npft),          ! maximum number of days past planting allowed before auto-shift to grain fill
!            mxmat(npft),           ! maximum number of days allowed past planting for physiological maturity to be reached
!            fleafi(npft),          ! initial fraction of aboveground allocation going to leaf before grain fill begins
!            fleaf(npft),           ! fraction of aboveground allocation going to leaf before grain fill begins
!            cgrain(npft),          ! fraction of grain dry matter that is carbon
!            convfact(npft),        ! factor converting kg/m2 of C in grain to bu/acre value
!            maxhi(npft),           ! maximum harvest index
!            fyield(npft),          ! fraction of C allocated to grain pool that is actually seed
!            fnlfmx(npft),          ! maximum amount of N allowed in leaf at end of growing season
!            fngrmx(npft),          ! maximum amount of N allowed in grain at end of growing season
!            sratio(npft),          ! leaf:stem N allocation ratio
!            rratio(npft),          ! leaf:root N allocation ratio
!            fnopt(npft),           ! minimum leaf nitrogen content that doesn't experience N stress
!            bfact(npft)            ! coefficient in LAI curve
!            grainmoisture(npft)    ! grain moisture, must be between 0.08 and 0.16
!
      real*8, dimension(:), allocatable :: grnwht, fleafiwht, fleafwht, mgddgf,        &
                                           mxmatwht, fnlfmxw, fngrmxw, fnoptw,         &
                                           cfertmaize, cfertsoy, cfertsgc, cfertwheat
!
!            grnwht(2),             ! fraction of annual gdd (to reach phys. mat.) for wheat leaf emergence
!            fleafiwht(2),          ! initial fraction of aboveground allocation going to leaf before grain fill begins
!            fleafwht(2),           ! fraction of aboveground allocation going to leaf before grain fill begins
!            mgddgf(2),             ! maximum gdd allowed past grain fill inititiation for phys. maturity
!            mxmatwht(2),           ! maximum number of days allowed past planting for physiological maturity to be reached
!            fnlfmxw(2),            ! leaf nitrogen maximum allowed for wheat
!            fngrmxw(2),            ! grain nitrogen maximum allowed for wheat
!            fnoptw(2)              ! optimum leaf nitrogen fraction for wheat
!            cfertmaize(2000),&     ! estimated  annual average N-fertilizer applied to maize   crops across US Mississippi basin (1950-00)
!            cfertsoy(2000), &      ! estimated  annual average N-fertilizer applied to soybean crops across US Mississippi Basin (1950-00)
!            cfertsgc(2000), &      ! estimated  annual average N-fertilizer applied to sugarcane crops across US Mississippi basin (1950-00)
!            cfertwheat(2000),      ! estimated  annual average N-fertilizer applied to wheat   crops across US Mississippi Basin (1950-00)
!
end module inland_comcrop
