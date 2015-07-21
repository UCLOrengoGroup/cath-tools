/*************************************************************************

   Program:    
   File:       OriginPDB.c
   
   Version:    V1.1R
   Date:       22.02.94
   Function:   Move a PDB linked list to the origin
   
   Copyright:  (c) SciTech Software 1993-4
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
   V1.1  22.02.94 Changed NULL check to any coordinate not 9999.0

*************************************************************************/
/* Includes
*/
#include <math.h>

#include "MathType.h"
#include "pdb.h"
#include "macros.h"

/************************************************************************/
/* Defines
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/* Variables global to this file only
*/

/***************************************************************************/
/*>void OriginPDB(PDB *pdb)
   ------------------------
   I/O:   PDB  *pdb    PDB linked list to move

   Moves a PDB linked list to the origin, ignoring NULL coordinates.

   01.10.92 Original
   22.02.94 Changed NULL check to any coordinate not 9999.0
   11.03.94 Changed NULL check to >9998.0 Added cast to REAL
*/
void OriginPDB(PDB *pdb)
{
   PDB   *p;
   VEC3F cg;

   GetCofGPDB(pdb,&cg);

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->x < (REAL)9999.0 || p->y < (REAL)9999.0 || p->z < (REAL)9999.0)
      {
         p->x -= cg.x;
         p->y -= cg.y;
         p->z -= cg.z;
      }
   }
}

