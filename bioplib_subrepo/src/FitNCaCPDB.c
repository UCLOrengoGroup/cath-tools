/*************************************************************************

   Program:    
   File:       FitNCaCPDB.c
   
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
/*>BOOL FitNCaCPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
   ----------------------------------------------------------
   Input:   PDB  *ref_pdb     Reference PDB linked list
   I/O:     PDB  *fit_pdb     Mobile PDB linked list
   Output:  REAL rm[3][3]     Rotation matrix (May be input as NULL).
   Returns: BOOL              Success

   Fits two PDB linked lists using only the CA atoms. 

   Actually fits fit_pdb onto ref_pdb and also returns the rotation 
   matrix. This may be NULL if these data are not required.

   12.12.01 Original based on FitCaPDB()   By: ACRM
*/
BOOL FitNCaCPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
{
   REAL  RotMat[3][3];
   COOR  *ref_coor   = NULL,
         *fit_coor   = NULL;
   VEC3F ref_bb_CofG,
         fit_bb_CofG,
         tvect;
   int   NCoor       = 0,
         i, j,
         natoms;
   BOOL  RetVal;
   PDB   *ref_bb_pdb = NULL,
         *fit_bb_pdb = NULL;
   char  *sel[4];

   /* First extract only the backbone (BB) atoms                         */
   SELECT(sel[0], "N   ");
   SELECT(sel[1], "CA  ");
   SELECT(sel[2], "C   ");
   if((sel[0]==NULL)||(sel[1]==NULL)||(sel[2]==NULL)||(sel[3]==NULL))
      return(FALSE);
   if((ref_bb_pdb = SelectAtomsPDB(ref_pdb, 3, sel, &natoms))==NULL)
      RetVal = FALSE;
   if((fit_bb_pdb = SelectAtomsPDB(fit_pdb, 3, sel, &natoms))==NULL)
      RetVal = FALSE;
   free(sel[0]);
   free(sel[1]);
   free(sel[2]);
   
   /* If we succeeded in building our BB PDB linked lists...            */
   if(RetVal)
   {
      /* Get the CofG of the BB structures and the original mobile      */
      GetCofGPDB(ref_bb_pdb, &ref_bb_CofG);
      GetCofGPDB(fit_bb_pdb, &fit_bb_CofG);
      
      /* Move them both to the origin                                   */
      OriginPDB(ref_bb_pdb);
      OriginPDB(fit_bb_pdb);
      
      /* Create coordinate arrays, checking numbers match               */
      NCoor = GetPDBCoor(ref_bb_pdb, &ref_coor);
      if(GetPDBCoor(fit_bb_pdb, &fit_coor) != NCoor)
      {
         RetVal = FALSE;
      }
      else
      {
         /* Can't fit with fewer than 3 coordinates                     */
         if(NCoor < 3)
         {
            RetVal = FALSE;
         }
         else
         {
            /* Everything OK, go ahead with the fitting                 */
            if(!matfit(ref_coor,fit_coor,RotMat,NCoor,NULL,FALSE))
            {
               RetVal = FALSE;
            }
            else
            {
               /* Apply the operations to the true coordinates          */
               tvect.x = (-fit_bb_CofG.x);
               tvect.y = (-fit_bb_CofG.y);
               tvect.z = (-fit_bb_CofG.z);
               TranslatePDB(fit_pdb, tvect);
               ApplyMatrixPDB(fit_pdb, RotMat);
               TranslatePDB(fit_pdb, ref_bb_CofG);
            }
         }
      }
   }
   
   /* Free the coordinate arrays and BB PDB linked lists                */
   if(ref_coor)   free(ref_coor);
   if(fit_coor)   free(fit_coor);
   if(ref_bb_pdb) FREELIST(ref_bb_pdb, PDB);
   if(fit_bb_pdb) FREELIST(fit_bb_pdb, PDB);
         
   /* Fill in the rotation matrix for output, if required               */
   if(RetVal && (rm!=NULL))
   {
      for(i=0; i<3; i++)
         for(j=0; j<3; j++)
            rm[i][j] = RotMat[i][j];
   }

   return(RetVal);
}

