/*************************************************************************

   Program:    
   File:       GetCofGPDBSCRange.c
   
   Version:    V1.1R
   Date:       03.10.94
   Function:   Find CofG of a PDB linked list
   
   Copyright:  (c) SciTech Software 1992-4
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
   V1.0  01.10.92 Original
   V1.1  03.10.94 Added GetCofGPDBRange(), FindCofGPDBSCRange() and 
                  fixed NULL coord search in GetCofGPDB()

*************************************************************************/
/* Includes
*/
#include <math.h>

#include "pdb.h"
#include "macros.h"
#include "MathType.h"

/************************************************************************/
/* Defines
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/* Variables global to this file only
*/


/************************************************************************/
/*>void GetCofGPDBSCRange(PDB *start, PDB *stop, VEC3F *cg)
   --------------------------------------------------------
   Input:    PDB   *start     Start of region of interest in PDB list
             PDB   *stop      Beginning of next residue
   Output:   VEC3F *cg        Centre of geometry of specified region

   Find CofG of a range in a PDB linked list, ignoring NULL coordinates
   (specified as all coords==9999.000) and backbone (N,CA,C,O).
   For Glycine, returns the CA coordinates.

   03.10.94 Original    By: ACRM
*/
void GetCofGPDBSCRange(PDB *start, PDB *stop, VEC3F *cg)
{
   int natom;
   PDB *p, *ca;

   cg->x = 0.0;   
   cg->y = 0.0;   
   cg->z = 0.0;   
   natom = 0;

   for(p=start; p!=NULL && p!=stop; NEXT(p))
   {
      if(!strncmp(p->atnam,"CA  ",4))
         ca = p;
      
      if(p->x < 9999.0 || p->y < 9999.0 || p->z < 9999.0)
      {
         if(strncmp(p->atnam,"N   ",4) &&
            strncmp(p->atnam,"CA  ",4) &&
            strncmp(p->atnam,"C   ",4) &&
            strncmp(p->atnam,"O   ",4))
         {
            cg->x += p->x;
            cg->y += p->y;
            cg->z += p->z;
            natom++;
         }
      }
   }

   if(natom)
   {
      cg->x /= natom;
      cg->y /= natom;
      cg->z /= natom;
   }
   else
   {
      cg->x = ca->x;
      cg->y = ca->y;
      cg->z = ca->z;
   }
}

