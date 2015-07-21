/*************************************************************************

   Program:    
   File:       ApMatPDB.c
   
   Version:    V1.0R
   Date:       August 1993
   Function:   
   
   Copyright:  (c) SciTech Software 1993
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

*************************************************************************/
/* Includes
*/
#include <math.h>
#include "MathType.h"
#include "pdb.h"
#include "matrix.h"
#include "macros.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/* Variables global to this file only
*/

/************************************************************************/
/*>void ApplyMatrixPDB(PDB *pdb, REAL matrix[3][3])
   ------------------------------------------------
   I/O:   PDB  *pdb          PDB linked list
   Input: REAL matrix[3][3]  Matrix to apply

   Apply a rotation matrix to a PDB linked list.

   22.07.93 Original (old RotatePDB())   By: ACRM
*/
void ApplyMatrixPDB(PDB  *pdb,
                    REAL matrix[3][3])
{
   PDB   *p;
   VEC3F incoords,
         outcoords;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->x != 9999.0 && p->y != 9999.0 && p->z != 9999.0)
      {
         incoords.x = p->x;
         incoords.y = p->y;
         incoords.z = p->z;
         MatMult3_33(incoords,matrix,&outcoords);
         p->x = outcoords.x;
         p->y = outcoords.y;
         p->z = outcoords.z;
      }
   }
}

