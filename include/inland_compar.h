/*
 * Compile-time main parameter definitions for the INLAND model.
 */

#ifndef __INCLUDED_COMPAR_H
#define __INCLUDED_COMPAR_H

#include "inland_config.h"

#ifndef SINGLE_POINT_MODEL
/* Definitions for the Grid model (2D) */
#define PLON 360 /* nlon */
#define PLAT 180 /* nlat */
/* LPT (npoi) not used anymore */
#define XRES 1.0 /* xres */
#define YRES 1.0 /* yres */

#else /* SINGLE_POINT_MODEL */
/* Default definitions for the Single Point model (0D) */
#define PLON 2 /* nlon */
#define PLAT 2 /* nlat */
/* LPT (npoi) not used anymore */
#define XRES 0.083333 /* xres */
#define YRES 0.083333 /* yres */

#endif /* SINGLE_POINT_MODEL */

/* General definitions */
#define NUMBANDS 2 /* nband */
#define NUMSOILLAYERS 11 /* nsoilay */
#define NUMSNOWLAYERS 3 /* nsnolay */
#ifndef SINGLE_POINT_MODEL
#define NUMPFT 16 /* npft */
#else
#define NUMPFT 12 /* npft */
#endif /* SINGLE_POINT_MODEL */
#define SCPFT 13 /* scpft */
#define ECPFT 16 /* ecpft */

#define PISTOS 3.1415927 /* PI */
#endif /* ndef __INCLUDED_COMPAR_H */
