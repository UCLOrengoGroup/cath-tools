/*************************************************************************

   Program:    
   File:       RotPDB.c
   
   Version:    V1.0R
   Date:       August 1993
   Function:   Rotate a PDB linked list
   
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
   Rotate a PDB linked list. Moves the structure to the origin first,
   applies the rotation and moves back from the origin.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdlib.h>

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


/************************************************************************/
/*>void RotatePDB(PDB *pdb, REAL matrix[3][3])
   -------------------------------------------
   I/O:    PDB   *pdb          PDB linked list to rotate
   Input:  REAL  matrix[3][3]  Rotation matrix

   Rotates a PDB linked list using ApplyMatrixPDB() which ignores 
   coordinates of 9999.0. The structure is moved to the origin, the 
   matrix is applied and the structure is moved back.

   30.09.92 Original
   01.10.92 Added check on NULL coordinates
   22.07.93 Moves to origin first; calls ApplyMatrixPDB() to do the work
*/
void RotatePDB(PDB *pdb, REAL matrix[3][3])
{
   VEC3F CofG;
         
   GetCofGPDB(pdb, &CofG);
   OriginPDB(pdb);
   
   ApplyMatrixPDB(pdb, matrix);
   
   TranslatePDB(pdb, CofG);
}

