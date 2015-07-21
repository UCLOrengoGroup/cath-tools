/*************************************************************************

   Program:    
   File:       WriteCrystPDB.c
   
   Version:    V1.0
   Date:       12.10.95
   Function:   
   
   Copyright:  (c) SciTech Software 1993-1995
   Author:     Dr. Andrew C. R. Martin
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

**************************************************************************

   Revision History:
   =================
   V1.0R 12.10.05 Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include "MathType.h"

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
/*>void WriteCrystPDB(FILE *fp, VEC3F UnitCell, VEC3F CellAngles,
                      char *spacegroup,
                      REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4])
   -----------------------------------------------------------------
   Input:   FILE    *fp            Output file pointet
            VEC3F   UnitCell       The unit cell dimensions
            VEC3F   CellAngles     The unit cell angles
            char    *spacegroup    The crystal's space group
            REAL    OrigMatrix     The origin matrix
            REAL    ScaleMatrix    The scale matrix

   Write crystal parameters (unit cell, space group, origin and scale
   matrices) to a PDB file.

   12.10.95 Original    By: ACRM
   17.10.95 Corrected %lf to %f in fprintf()s
*/
void WriteCrystPDB(FILE *fp, VEC3F UnitCell, VEC3F CellAngles,
                   char *spacegroup,
                   REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4])
{
   int i;
   
   CellAngles.x *= 180.0/PI;
   CellAngles.y *= 180.0/PI;
   CellAngles.z *= 180.0/PI;

   fprintf(fp,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %14s\n",
           UnitCell.x, UnitCell.y, UnitCell.z,
           CellAngles.x, CellAngles.y, CellAngles.z, spacegroup);
   for(i=0; i<3; i++)
   {
      fprintf(fp,"ORIGX%1d    %10.6f%10.6f%10.6f%15.5f\n", i+1,
              OrigMatrix[i][0],
              OrigMatrix[i][1],
              OrigMatrix[i][2],
              OrigMatrix[i][3]);
   }
   for(i=0; i<3; i++)
   {
      fprintf(fp,"SCALE%1d    %10.6f%10.6f%10.6f%15.5f\n", i+1,
              ScaleMatrix[i][0],
              ScaleMatrix[i][1],
              ScaleMatrix[i][2],
              ScaleMatrix[i][3]);
   }
}

