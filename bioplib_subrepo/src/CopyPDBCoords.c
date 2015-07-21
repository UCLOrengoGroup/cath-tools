/*************************************************************************

   Program:    
   File:       CopyPDBCoords.c
   
   Version:    V1.10R
   Date:       08.10.99
   Function:   PDB linked list manipulation
   
   Copyright:  (c) SciTech Software 1992-6
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      andrew@stagleys.demon.co.uk
               
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
   V1.0  22.02.94 Original release
   V1.1  23.05.94 Added FindNextChainPDB()
   V1.2  05.10.94 KillSidechain() uses BOOL rather than int
   V1.3  24.07.95 Added TermPDB()
   V1.4  25.07.95 Added GetPDBChainLabels()
   V1.5  26.09.95 Fixed bug in TermPDB()
   V1.6  12.10.95 Added DupePDB(), CopyPDBCoords()
   V1.7  23.10.95 Moved FindResidueSpec() to ParseRes.c
   V1.8  10.01.96 Added ExtractZonePDB()
   V1.9  14.03.96 Added FindAtomInRes()
   V1.10 08.10.99 Initialised some variables

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdlib.h>

#include "MathType.h"
#include "SysDefs.h"
#include "pdb.h"
#include "macros.h"
#include "general.h"

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
/*>BOOL CopyPDBCoords(PDB *out, PDB *in)
   -------------------------------------
   Input:   PDB  *in      Input PDB linked list
   Output:  PDB  *out     Output PDB linked list
   Returns: BOOL          Success?

   Applies the coordinates of `in' to `out'. Assumes that the structures
   are equivalent with identical atom ordering. Makes a simple check on
   resnam and atnam at each position.

   11.10.95 Original   By: ACRM
*/
BOOL CopyPDBCoords(PDB *out, PDB *in)
{
   PDB *p, *q;
   
   for(p=in, q=out; p!=NULL && q!=NULL; NEXT(p), NEXT(q))
   {
      if(strncmp(p->atnam,  q->atnam,  4) ||
         strncmp(p->resnam, q->resnam, 4))
         return(FALSE);
      
      q->x = p->x;
      q->y = p->y;
      q->z = p->z;
   }

   if(p!=NULL || q!=NULL)
      return(FALSE);
   
   return(TRUE);
}


