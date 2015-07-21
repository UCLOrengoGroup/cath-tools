/*************************************************************************

   Program:    
   File:       CalcCterCoords.c
   
   Version:    V1.3R
   Date:       13.11.96
   Function:   Calcualtes C-terminal oxygen coordinates.
   
   Copyright:  (c) SciTech Software 1994-6
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      andrew@bioinf.org.uk
               
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

   BOOL CalcCterCoords(PDB *p, PDB *ca_p, PDB *c_p, PDB *o_p)
   ----------------------------------------------------------
   Calculates the coordinates for a second oxygen (p) given the
   3 antecedent atoms. Normally called from FixCterPDB()

**************************************************************************

   Revision History:
   =================
   V1.0  24.08.94 Original    By: ACRM
   V1.1  05.10.94 Removed unused variables
   V1.2  12.11.96 If any of the antecedant coordinates are undefined, set
                  the terminal oxygen to NULL coordinates
   V1.3  13.11.96 Also checks for missing CA,C and O1 records

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SysDefs.h"
#include "MathType.h"
#include "pdb.h"
#include "macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF      160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>BOOL CalcCterCoords(PDB *p, PDB *ca_p, PDB *c_p, PDB *o_p)
   ----------------------------------------------------------
   I/O:     PDB  *p     OT2 PDB record whose coords are to be fixed
   Input:   PDB  *ca_p  Antecedent CA PDB pointer
            PDB  *c_p   Antecedent C PDB pointer
            PDB  *o_p   Antecedent O PDB pointer
   Returns: BOOL        Success

   This routine actually calculates the CTER OT2 coords.

   15.07.90 Original
*/
BOOL CalcCterCoords(PDB *p, PDB *ca_p, PDB *c_p, PDB *o_p)
{
   REAL  gr    = 1.3,
         alpha = 120.0*PI/180.0,
         cosa,  sina,  scalpr,
         x21,   y21,   z21,   d21,
         x23,   y23,   z23,   d23,
         x32,   y32,   z32,
         xp23,  yp23,  zp23,  rp23,
         xh,    yh,    zh,
         xp,    yp,    zp,
         xv,    yv,    zv;
         
   if(ca_p==NULL || c_p==NULL || o_p==NULL)
   {
      return(FALSE);
   }
   
   x21 = c_p->x - ca_p->x;
   y21 = c_p->y - ca_p->y;
   z21 = c_p->z - ca_p->z;
   d21 = (REAL)sqrt((double)(x21*x21 + y21*y21 + z21*z21));
   
   x23 = c_p->x - o_p->x;
   y23 = c_p->y - o_p->y;
   z23 = c_p->z - o_p->z;
   d23 = (REAL)sqrt((double)(x23*x23 + y23*y23 + z23*z23));

   cosa = (REAL)cos((double)alpha);
   sina = (REAL)sin((double)alpha);

   x32 = o_p->x - c_p->x;
   y32 = o_p->y - c_p->y;
   z32 = o_p->z - c_p->z;

   scalpr = (x21*x32 + y21*y32 + z21*z32)/d21;
   xh = x21 / d21;
   yh = y21 / d21;
   zh = z21 / d21;

   xp = scalpr * xh;
   yp = scalpr * yh;
   zp = scalpr * zh;

   xp23 = xp + x23;
   yp23 = yp + y23;
   zp23 = zp + z23;

   rp23 = (REAL)sqrt((double)(xp23*xp23 + yp23*yp23 + zp23*zp23));
   xv = xp23/rp23;
   yv = yp23/rp23;
   zv = zp23/rp23;

   p->x = c_p->x + gr*(-cosa*xh + sina*xv);
   p->y = c_p->y + gr*(-cosa*yh + sina*yv);
   p->z = c_p->z + gr*(-cosa*zh + sina*zv);
   
   return(TRUE);
}

