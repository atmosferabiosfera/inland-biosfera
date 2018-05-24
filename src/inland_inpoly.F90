#include "inland_config.h"
! ---------------------------------------------------------------------------
     logical function inpoly (poly, npoints, xt, yt)
!
! adapted by cjk 01/11/01
! ---------------------------------------------------------------------------
!
!                            INPOLY
! Function to tell if a point is inside a polygon or not.
! --------------------------------------------------------------------------
! Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.
!
! Please feel free to use this source code for any purpose, commercial
! or otherwise, as long as you don't restrict anyone else's use of
! this source code.  Please give credit where credit is due.
!
! Point-in-polygon algorithm, created especially for World-Wide Web
! servers to process image maps with mouse-clickable regions.
!
! Home for this file:  http://www.gcomm.com/develop/inpoly.c
!
!                 6/19/95 - Bob Stein & Craig Yap
!                       stein@gcomm.com
!                      craig@cse.fau.edu
! --------------------------------------------------------------------------
! Modified by:
! Aris Gerakis, apr. 1998: 1.  translated to Fortran
!                          2.  made it work with real coordinates
!                          3.  now resolves the case where point falls
!                              on polygon border.
! Aris Gerakis, nov. 1998: Fixed error caused by hardware arithmetic
! Aris Gerakis, july 1999: Now all borderline cases are valid
! --------------------------------------------------------------------------
! Glossary:
! function inpoly: true=inside, false=outside (is target point inside
!                    a 2D polygon?)
!   poly(*,2):  polygon points, [0]=x, [1]=y
!   npoints: number of points in polygon
!   xt: x (horizontal) of target point
!   yt: y (vertical) of target point
! --------------------------------------------------------------------------
!
! input variables
!
     integer ::  npoints, &
                 xt,      &
                 yt
!
     real*8 ::   poly(7, 2)
!
! local variables
!
     real*8 ::  xnew,   &
                ynew,   &
                xold,   &
                yold,   &
                x1,     &
                y1,     &
                x2,     &
                y2
!
     integer ::  i
!
     logical :: inside, on_border

     inside = .false.
     on_border = .false.
!
     if (npoints .lt. 3)  then
        inpoly = .false.
        return
     end if
!
     xold = poly(npoints,1)
     yold = poly(npoints,2)

     do 300  i = 1 , npoints
        xnew = poly(i,1)
        ynew = poly(i,2)

        if (xnew .gt. xold)  then
           x1 = xold
           x2 = xnew
           y1 = yold
           y2 = ynew
        else
           x1 = xnew
           x2 = xold
           y1 = ynew
           y2 = yold
        end if

! the outer IF is the 'straddle' test and the 'vertical border' test.
! the inner IF is the 'non-vertical border' test and the 'north' test.

! the first statement checks whether a north pointing vector crosses
! (stradles) the straight segment.  There are two possibilities, depe-
! nding on whether xnew < xold or xnew > xold.  The '<' is because edge
! must be "open" at left, which is necessary to keep correct count when
! vector 'licks' a vertix of a polygon.

        if ((xnew .lt. xt .and. xt .le. xold) .or. &
            (.not. xnew .lt. xt .and. .not. xt .le. xold)) then
!
! the test point lies on a non-vertical border:
!
           if ((yt-y1)*(x2-x1) .eq. (y2-y1)*(xt-x1)) then
              on_border = .true.
!
! check if segment is north of test point.  If yes, reverse the
! value of INSIDE.  The +0.001 was necessary to avoid errors due
! arithmetic (e.g., when clay = 98.87 and sand = 1.13):
!
           elseif ((yt-y1)*(x2-x1) .lt. (y2-y1)*(xt-x1) + 0.001) then
              inside = .not.inside ! cross a segment
           endif
!
! this is the rare case when test point falls on vertical border or
! left edge of non-vertical border. The left x-coordinate must be
! common.  The slope requirement must be met, but also point must be
! between the lower and upper y-coordinate of border segment.  There
! are two possibilities,  depending on whether ynew < yold or ynew >
! yold:
!
        elseif ((xnew .eq. xt .or. xold .eq. xt) .and. (yt-y1)*(x2-x1) &
                 .eq. (y2-y1)*(xt-x1) .and. ((ynew .le. yt .and. yt    &
                 .le. yold) .or. (.not. ynew .lt. yt .and. .not. yt    &
                 .lt. yold))) then
                    on_border = .true.
        endif
!
        xold = xnew
        yold = ynew
!
 300    continue  
!
! If test point is not on a border, the function result is the last state
! of INSIDE variable.  Otherwise, INSIDE doesn't matter.  The point is
! inside the polygon if it falls on any of its borders:
!
     if (.not. on_border) then
        inpoly = inside
     else
        inpoly = .true.
     endif
!
      end function inpoly
