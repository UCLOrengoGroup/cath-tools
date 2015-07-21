/*************************************************************************

   Program:    
   File:       angle.c
   
   Version:    V1.5R
   Date:       27.03.95
   Function:   Calculate angles, torsions, etc.
   
   Copyright:  (c) SciTech Software 1993-5
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
   These routines return angles and torsion angles. The definition of a
   torsion angle is the chemical definition:
   i.e. Assuming the atoms are co-planar:
            A---B          A---B
                | = 0.0        |  = 180.0
                |              |
            D---C              C---D

**************************************************************************

   Usage:
   ======
   REAL angle(xi,yi,zi,xj,yj,zj,xk,yk,zk)
   Input:   REAL     xi,yi,zi    Input coordinates
                     xj,yj,zj
                     xk,yk,zk
   Returns: REAL                 The angle between the 3 atoms


   REAL phi(xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl)
   Input:   REAL     xi,yi,zi    Input coordinates
                     xj,yj,zj
                     xk,yk,zk
                     xl,yl,zl
   Returns: REAL                 The torsion angle between the 4 atoms


   REAL simpleang(ang)
   Input:   REAL     ang         An angle
   Returns: REAL                 Simplified angle
   
   Simplifies a signed angle to an unsigned angle <=2*PI


   REAL TrueAngle(REAL opp, REAL adj)
   Input:   REAL     opp         Length of opposite side
            REAL     adj         Length of adjacent side
   Returns: REAL                 The angle from 0 to 2PI

   Returns the true positive angle between 0 and 2PI given the opp and
   adj lengths

**************************************************************************

   Revision History:
   =================
   V1.0  07.02.91 Original
   V1.1  17.02.91 Corrected comments to new standard and added phi()
   V1.2  04.03.91 angle() and phi() now return _correct_ values!
   V1.3  01.06.92 ANSIed
   V1.4  08.12.92 Changed abs() to ABS() from macros.h
   V1.5  27.03.95 Added TrueAngle()

*************************************************************************/
/* Includes
*/
#include <math.h>

#include "MathType.h"
#include "macros.h"

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
/*>REAL angle(REAL xi,REAL yi,REAL zi,REAL xj,REAL yj,
              REAL zj,REAL xk,REAL yk,REAL zk)
   ---------------------------------------------------
   Input:   REAL    xi,yi,zi    Input coordinates
                    xj,yj,zj
                    xk,yk,zk
   Returns: REAL                The angle between the 3 atoms

   Calculates the angle between three sets of coordinates

   07.02.89 Original    By: ACRM
   04.03.91 Fixed return value
   16.06.93 Changed float to REAL
*/
REAL angle(REAL xi,
           REAL yi,
           REAL zi,
           REAL xj,
           REAL yj,
           REAL zj,
           REAL xk,
           REAL yk,
           REAL zk)
{
   REAL qx,qy,qz,sq,px,py,pz,sp,cosa2,a2;

   px = xi - xj;
   py = yi - yj;
   pz = zi - zj;
   sp = sqrt(px * px + py * py + pz * pz);

   qx = xk - xj;
   qy = yk - yj;
   qz = zk - zj;
   sq = sqrt(qx * qx + qy * qy + qz * qz);

   cosa2 = (qx * px + qy * py + qz * pz) / (sp * sq);
   a2 = acos(cosa2);

   return(a2);
}

