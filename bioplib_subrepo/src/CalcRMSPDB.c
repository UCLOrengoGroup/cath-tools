/*************************************************************************

   Program:    
   File:       CalcRMSPDB.c
   
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
/*>REAL CalcRMSPDB(PDB *pdb1, PDB *pdb2)
   -------------------------------------
   Input:   PDB  *pdb1   First PDB linked list
            PDB  *pdb2   Second PDB linked list
   Returns: REAL         RMS deviation

   Calculate the RMS deviation between two fitted PDB linked lists. The
   two lists must contain equivalent structures (same atom types in same
   order). No checks are made on this.

   11.03.94 Original    By: ACRM
*/
REAL CalcRMSPDB(PDB *pdb1, PDB *pdb2)
{
   PDB *p, *q;
   int count = 0;
   REAL dist = (REAL)0.0,
        rms;
   
   for(p=pdb1, q=pdb2; p!=NULL && q!=NULL; NEXT(p), NEXT(q))
   {
      dist += DISTSQ(p,q);
      count++;
   }

   rms = (REAL)((count)?sqrt((double)(dist/(REAL)count)):0.0);
   return(rms);
}


