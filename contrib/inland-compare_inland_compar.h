/*
 * Compile-time main parameter definitions for the INLAND model.
 */

#ifndef __INCLUDED_COMPAR_H
#define __INCLUDED_COMPAR_H

#include "inland_config.h"

#ifndef SINGLE_POINT_MODEL

/* Definitions for the Grid model (2D) */
/* longitude and latitude dimension of domain */
#define PLON 720 /* nlon */
#define PLAT 360 /* nlat */
#define LPT 1680 /* number of land points (npoi) */
#define XRES 0.5 /* xres */
#define YRES 0.5 /* yres */

#else /* SINGLE_POINT_MODEL */

/* Default definitions for the Single Point model (0D) */
#define PLON 2 /* nlon */
#define PLAT 2 /* nlat */
#define LPT 1 /* npoi */
#define XRES 0.083333 /* xres */
#define YRES 0.083333 /* yres */

#endif /* SINGLE_POINT_MODEL */

/* General definitions */
#define NUMBANDS 2 /* nband */
#define NUMSOILLAYERS 6 /* nsoilay */
#define NUMSNOWLAYERS 3 /* nsnolay */
#define NUMPFT 12 /* npft */

#define PISTOS 3.1415927 /* PI */
#endif /* ndef __INCLUDED_COMPAR_H */
