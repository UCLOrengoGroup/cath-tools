/*************************************************************************

   Program:    
   File:       angle.h
   
   Version:    V1.7R
   Date:       06.09.96
   Function:   Include file for angle functions
   
   Copyright:  (c) SciTech Software 1993-6
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      martin@biochem.ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.5  27.03.95 
   V1.6  09.07.96 Added TorToCoor()
   V1.7  06.09.96 Added includes

*************************************************************************/
/* Includes
*/
#include "MathType.h"
#include "SysDefs.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/

#ifndef _ANGLE_H
#define _ANGLE_H

REAL angle(REAL xi, REAL yi, REAL zi, REAL xj, REAL yj, REAL zj,
           REAL xk, REAL yk, REAL zk);
REAL phi(REAL xi, REAL yi, REAL zi, REAL xj, REAL yj, REAL zj,
         REAL xk, REAL yk, REAL zk, REAL xl, REAL yl, REAL zl);
REAL simpleangle(REAL ang);
REAL TrueAngle(REAL opp, REAL adj);
BOOL TorToCoor(VEC3F ant1, VEC3F ant2, VEC3F ant3, 
               REAL bond, REAL theta, REAL torsion,
               VEC3F *coords);

#endif
