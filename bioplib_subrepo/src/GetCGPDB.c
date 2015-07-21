/*************************************************************************

   Program:    
   File:       GetCGPDB.c
   
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
/*>void GetCofGPDB(PDB *pdb,VEC3F *cg)
   -----------------------------------
   Input:    PDB   *pdb       Start of PDB linked list
   Output:   VEC3F *cg        Centre of geometry of specified region

   Finds the CofG of a PDB linked list, ignoring NULL coordinates.

   01.10.92 Original
   03.10.94 Fixed NULL coordinate ignoring
*/
void GetCofGPDB(PDB   *pdb,
                VEC3F *cg)
{
   int natom;
   PDB *p;

   cg->x = 0.0;   
   cg->y = 0.0;   
   cg->z = 0.0;   
   natom = 0;
   for(p=pdb; p; NEXT(p))
   {
      if(p->x < 9999.0 || p->y < 9999.0 || p->z < 9999.0)
      {
         cg->x += p->x;
         cg->y += p->y;
         cg->z += p->z;
         natom++;
      }
   }
   cg->x /= natom;
   cg->y /= natom;
   cg->z /= natom;
}

