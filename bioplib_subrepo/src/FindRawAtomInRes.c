/*************************************************************************

   Program:    
   File:       FindRawAtomInRes.c
   
   Version:    V1.12R
   Date:       03.06.05
   Function:   PDB linked list manipulation
   
   Copyright:  (c) SciTech Software 1992-2005
   Author:     Dr. Andrew C. R. Martin
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
   V1.11 28.02.01 Added FindRawAtomInRes()
   V1.12 03.06.05 Compares 4 rather than 5 characters

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
/*>PDB *FindRawAtomInRes(PDB *pdb, char *atnam_in)
   -----------------------------------------------
   Input:   PDB    *pdb         The beginning of a residue in a PDB 
                                linked list
            char   *atnam_in    An atom name to search for (doesn't need
                                to be space-padded)
   Returns: PDB    *            Pointer to required atom, NULL if not
                                found

   28.02.01 Original based on FindAtomInRes()  By: ACRM
   03.06.05 Now compares 4 characters rather than 5
*/
PDB *FindRawAtomInRes(PDB *pdb, char *atnam_in)
{
   PDB *end,
       *p;

   char atnam[8];
   
   /* First copy the specified atom name and pad to 5 chars             */
   strcpy(atnam,atnam_in);
   padterm(atnam,5);
   
   /* Find the end of this residue                                      */
   end = FindNextResidue(pdb);
   
   /* Search for the required atom                                      */
   for(p=pdb; p!=end; NEXT(p))
   {
      if(!strncmp(p->atnam_raw,atnam,4))
         return(p);
   }
   
   return(NULL);
}
