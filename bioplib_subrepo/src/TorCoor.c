/*************************************************************************

   Program:    
   File:       TorCoor.c
   
   Version:    V1.0
   Date:       08.07.96
   Function:   Calculate cartesian coordinates from torsion angle
   
   Copyright:  (c) SciTech Software 1989-96
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      andrew@stagleys.demon.co.uk
               
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

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdio.h>

#include "MathType.h"
#include "SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
#define ETA (1.07e-7)

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>BOOL TorToCoor(VEC3F ant1, VEC3F ant2, VEC3F ant3, 
                  REAL bond, REAL theta, REAL torsion,
                  VEC3F *coords)
   ---------------------------------------------------
   Input:   VEC3F     ant1      First antecedent atom coordinates
            VEC3F     ant2      Second antecedent atom coordinates
            VEC3F     ant3      Third antecedent atom coordinates
            REAL      bond      Bond length from ant3 to new atom
            REAL      theta     Bond angle ant2-ant3-new
            REAL      torsion   Torsion angle ant1-ant2-ant3-new
   Output:  VEC3F     *coords   Coordinates of new atom
   Returns: BOOL                TRUE if distance between atoms 2 and 3
                                is < ETA (1.07e-7)

   Calculates cartesian coordinates for an atom given the coordinates of
   three antecedant atoms and the bond length, angle and torsion angle

   08.07.96 Original By: ACRM based on FORTRAN code adapted from Bob
            Bruccoleri's code from CONGEN.
*/
BOOL TorToCoor(VEC3F ant1, VEC3F ant2, VEC3F ant3, 
               REAL bond, REAL theta, REAL torsion,
               VEC3F *coords)
{
   REAL x1,  y1,  z1,
        x2,  y2,  z2,
        x3,  y3,  z3,
        xx4, yy4, zz4,
        x2o, y2o, z2o,
        CosTheta, SinTheta,
        CosTor,   SinTor,
        DistYZ1,  InvDistYZ1,
        DistXZ2,  DistSqXZ2,
        Dist2,    InvDist2,
        y1o,      z1o,
        xz2o,     xx1,
        BondSinTheta,
        InvDistXZ2;
   BOOL RetVal = FALSE;

   /* Generate pretend position of new atom if all other atoms were
      easily lined up.
   */
   SinTheta = sin(PI-theta);
   CosTheta = cos(PI-theta);
   SinTor = sin(torsion);
   CosTor = cos(torsion);
   BondSinTheta = bond * SinTheta;
   coords->x = bond * CosTheta;
   coords->y = BondSinTheta * CosTor;
   coords->z = BondSinTheta * SinTor;

   /* Translate ant1 and ant2 such that ant3 is at the origin.          */
   x3 = ant3.x;
   y3 = ant3.y;
   z3 = ant3.z;
   x1 = ant1.x - x3;
   y1 = ant1.y - y3;
   z1 = ant1.z - z3;
   x2 = ant2.x - x3;
   y2 = ant2.y - y3;
   z2 = ant2.z - z3;

   /* Rotate ant1 by rotation of ant2 to the origin.                    */
   DistSqXZ2 = x2*x2 + z2*z2;
   Dist2     = sqrt(DistSqXZ2+y2*y2);
   DistXZ2   = sqrt(DistSqXZ2);
   
   if (Dist2 < ETA)
   {
      RetVal   = TRUE;
      InvDist2 = 1.0/ETA;
   }
   else
   {
      InvDist2 = 1.0/Dist2;
   }

   if (DistXZ2 < ETA)
   {
      xx1 = x1;
      x2o = 1.0;
      z2o = 0.0;
   }
   else
   {
      InvDistXZ2 = 1.0/DistXZ2;
      x2o        = x2*InvDistXZ2;
      z2o        = z2*InvDistXZ2;
      xx1        = x1*x2o + z1*z2o;
      z1         = z1*x2o - x1*z2o;
   }

   xz2o = DistXZ2 * InvDist2;
   y2o  = y2 * InvDist2;
   x1   = (-xx1*xz2o) - y1*y2o;
   y1   = xx1*y2o - y1*xz2o;
   
   /* Rotate new atom by inverse of rotation which takes the transformed 
      ant1 atom to the xy plane by rotation about the x axis.
   */
   DistYZ1    = sqrt(y1*y1 + z1*z1);
   InvDistYZ1 = 1.0 / DistYZ1;
   y1o        = y1 * InvDistYZ1;
   z1o        = z1 * InvDistYZ1;
   yy4        = y1o*coords->y - z1o*coords->z;
   zz4        = y1o*coords->z + z1o*coords->y;

   /* Rotate new atom by inverse of ant2 rotation to -x axis.           */
   xx4       = y2o*yy4 - xz2o*coords->x;
   coords->y = (-xz2o*yy4) - y2o*coords->x;
   coords->x = x2o*xx4 - z2o*zz4;
   coords->z = z2o*xx4 + x2o*zz4;

   /* Move back since we translated such that ant3 was at the origin    */
   coords->x += x3;
   coords->y += y3;
   coords->z += z3;

   return(RetVal);
}

