/*************************************************************************

   Program:    
   File:       CalcTertaHCoords.c
   
   Version:    V1.3R
   Date:       13.11.96
   Function:   Routines to add N-terminal hydrogens and C-terminal
               oxygens.
   
   Copyright:  (c) SciTech Software 1994-6
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      martin@biochem.ucl.ac.uk
               andrew@stagleys.demon.co.uk
               
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

   int CalcTetraHCoords(PDB *nter, COOR *coor)
   -------------------------------------------
   Calculate the coordinates for 3 tetrahedral hydrogens given a pointer
   to the residue onto which they are to be added. Normally called from
   AddNTerHs()

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
/*>int CalcTetraHCoords(PDB *nter, COOR *coor)
   -------------------------------------------
   Input:   PDB    *nter     Pointer to the N-terminus
   Output:  COOR   *coor     Array of hydrogen coordinates
   Returns: int              Number of hydrogens calculated (3)
                             0 if antecedant atoms missing

   Calculates coordinates for the extra hydrogens.

   23.08.94 Original    By: ACRM
*/
int CalcTetraHCoords(PDB *nter, COOR *coor)
{
   PDB *N  = NULL,
       *CA = NULL,
       *C  = NULL,
       *p,
       *end;
   REAL x21,     y21,     z21,     r21,
        x32,     y32,     z32,     r32,
        xh,      yh,      zh,      scalpr,
        xp,      yp,      zp,
        x21p,    y21p,    z21p,    r21p,
        xv,      yv,      zv,
        xs,      ys,      zs,
        cosax,   cosay,   cosaz,
        cosa,    sina,
        BondLen, sFac;

   cosa    = (REAL)cos((double)(70.6 * PI / 180.0));
   sina    = (REAL)sin((double)(70.6 * PI / 180.0));
   BondLen = (REAL)1.08;
   sFac    = (REAL)(sqrt((double)3.0) * 0.5);

   end = FindEndPDB(nter);

   /* Search for the antecedant atom pointers                           */
   for(p=nter; p!= end; NEXT(p))
   {
      if(!strncmp(p->atnam,"N   ",4) || !strncmp(p->atnam,"NT  ",4))
         N = p;
      else if(!strncmp(p->atnam,"CA  ",4))
         CA = p;
      else if(!strncmp(p->atnam,"C   ",4))
         C = p;
   }
   
   /* Check all were found                                              */
   if(N==NULL || CA==NULL || C==NULL)
   {
      fprintf(stderr,"Atom N,CA or C missing from residue: %s %c%d%c\n",
              p->resnam, p->chain[0], p->resnum, p->insert[0]);
      return(0);
   }

   /* Calculate coordinates for the 3 hydrogens                         */
   x21 = CA->x - C->x;
   y21 = CA->y - C->y;
   z21 = CA->z - C->z;
   r21 = (REAL)sqrt((double)(x21*x21 + y21*y21 + z21*z21));

   x32 = N->x - CA->x;
   y32 = N->y - CA->y;
   z32 = N->z - CA->z;
   r32 = (REAL)sqrt((double)(x32*x32 + y32*y32 + z32*z32));

   xh = x32/r32;
   yh = y32/r32;
   zh = z32/r32;

   scalpr = (x21*x32 + y21*y32 + z21*z32)/r32;

   xp = scalpr*xh;
   yp = scalpr*yh;
   zp = scalpr*zh;

   x21p = x21-xp;
   y21p = y21-yp;
   z21p = z21-zp;
   r21p = (REAL)sqrt((double)(x21p*x21p + y21p*y21p + z21p*z21p));

   xv = x21p/r21p;
   yv = y21p/r21p;
   zv = z21p/r21p;

   xs = yh*zv - zh*yv;
   ys = zh*xv - xh*zv;
   zs = xh*yv - yh*xv;

   cosax = cosa*xh;
   cosay = cosa*yh;
   cosaz = cosa*zh;
   
   coor[0].x = N->x + BondLen*(cosax + sina*xv);
   coor[0].y = N->y + BondLen*(cosay + sina*yv);
   coor[0].z = N->z + BondLen*(cosaz + sina*zv);
   
   coor[1].x = N->x + BondLen*(cosax + sina*(sFac*xs - 0.5*xv));
   coor[1].y = N->y + BondLen*(cosay + sina*(sFac*ys - 0.5*yv));
   coor[1].z = N->z + BondLen*(cosaz + sina*(sFac*zs - 0.5*zv));

   coor[2].x = N->x + BondLen*(cosax + sina*(-sFac*xs - 0.5*xv));
   coor[2].y = N->y + BondLen*(cosay + sina*(-sFac*ys - 0.5*yv));
   coor[2].z = N->z + BondLen*(cosaz + sina*(-sFac*zs - 0.5*zv));

   return(3);
}

