/*************************************************************************

   Program:    
   File:       FitPDB.c
   
   Version:    V1.3R
   Date:       14.03.96
   Function:   Fit two PDB linked lists. Also a weighted fit and support
               routines
   
   Copyright:  (c) SciTech Software 1993-6
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
   V1.0  01.03.94 Original release
   V1.1  11.03.94 Fixed bug in calls to matfit(). Had not been changed 
                  to reflect modification in MatMult3_33().
   V1.2  14.03.94 Fixed FitPDB(); wasn't filling in the output matrix
   V1.3  14.03.96 Added FitCaPDB()
                  Changed FitPDB() and FitCaCbPDB() to use 
                  ApplyMatrixPDB() rather than RotatePDB() since the PDB
                  linked lists are already at the origin

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "MathType.h"
#include "SysDefs.h"
#include "macros.h"
#include "fit.h"
#include "pdb.h"

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
/*>BOOL FitPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
   ------------------------------------------------------
   Input:   PDB  *ref_pdb     Reference PDB linked list
   I/O:     PDB  *fit_pdb     Mobile PDB linked list
   Output:  REAL rm[3][3]     Rotation matrix (May be input as NULL).
   Returns: BOOL              Success

   Fits two PDB linked lists. Actually fits fit_pdb onto ref_pdb and also
   returns the rotation matrix. This may be NULL if these data are not
   required.

   17.06.93 Original based on code from ProFit  By: ACRM
   11.03.94 Changed call to matfit(). Corrected to normal matrix.
   14.03.94 Actually fills in the rotation matrix (!). Restores original
            data if fitting failed.
   14.03.96 Changed to use ApplyMatrixPDB() rather than RotatePDB() since
            we are already at the origin
*/
BOOL FitPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
{
   REAL  RotMat[3][3];
   COOR  *ref_coor   = NULL,
         *fit_coor   = NULL;
   VEC3F ref_CofG,
         fit_CofG;
   int   NCoor       = 0,
         i, j;
   BOOL  RetVal;

   /* Get the CofG of the reference structure                           */
   GetCofGPDB(ref_pdb, &ref_CofG);
   GetCofGPDB(fit_pdb, &fit_CofG);

   /* Move them both to the origin                                      */
   OriginPDB(ref_pdb);
   OriginPDB(fit_pdb);
   
   /* Create coordinate arrays                                          */
   NCoor = GetPDBCoor(ref_pdb, &ref_coor);
   if(GetPDBCoor(fit_pdb, &fit_coor) != NCoor) 
   {
      /* Free and return if arrays don't match                          */
      if(ref_coor) free(ref_coor);
      if(fit_coor) free(fit_coor);
      return(FALSE);
   }
   
   /* Can't fit with fewer than 3 coordinates                           */
   if(NCoor < 3)
   {
      if(ref_coor) free(ref_coor);
      if(fit_coor) free(fit_coor);
      return(FALSE);
   }
   
   /* Everything OK, go ahead with the fitting                          */
   RetVal = matfit(ref_coor,fit_coor,RotMat,NCoor,NULL,FALSE);
   
   /* Now we can rotate the rotation list                               */
   if(RetVal)
   {
      ApplyMatrixPDB(fit_pdb, RotMat);
      TranslatePDB(fit_pdb, ref_CofG);
      TranslatePDB(ref_pdb, ref_CofG);
   }
   else
   {
      TranslatePDB(fit_pdb, fit_CofG);
      TranslatePDB(ref_pdb, ref_CofG);
   }

   /* Free the coordinate arrays                                        */
   free(ref_coor);
   free(fit_coor);

   /* Fill in the rotation matrix for output, if required               */
   if(rm!=NULL)
   {
      for(i=0; i<3; i++)
         for(j=0; j<3; j++)
            rm[i][j] = RotMat[i][j];
   }

   return(RetVal);
}

