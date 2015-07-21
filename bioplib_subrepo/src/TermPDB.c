/*************************************************************************

   Program:    
   File:       TermPDB.c
   
   Version:    V1.10
   Date:       08.10.99
   Function:   
   
   Copyright:  (c) SciTech Software 1992-6
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
#include "pdb.h"
#include "macros.h"

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
/*>PDB *TermPDB(PDB *pdb, int length)
   ----------------------------------
   Input:   PDB   *pdb         PDB linked list
            int   length       Number of residues after which to terminate
   Returns: PDB   *            Pointer to next residue after terminated
                               list. NULL if not enough residues in linked
                               list.

   Terminate a PDB linked list after length residues, returning a pointer
   to the next residue. 

   Note that the number of residues may cross chain boundaries.

   06.07.95 Original    By: ACRM
   26.09.95 Corrected update of resnum etc to use p-> not pdb-> (!!)
*/
PDB *TermPDB(PDB *pdb, int length)
{
   int  resnum,
        count;
   char insert,
        chain;
   PDB  *p,
        *prev = NULL;
   
   resnum = pdb->resnum;
   insert = pdb->insert[0];
   chain  = pdb->chain[0];

   for(p=pdb, count=1; p!=NULL; NEXT(p))
   {
      if((p->resnum    != resnum) ||
         (p->chain[0]  != chain)  ||
         (p->insert[0] != insert))
      {
         if(++count > length)
         {
            prev->next = NULL;
            return(p);
         }

         resnum = p->resnum;
         insert = p->insert[0];
         chain  = p->chain[0];
      }
      prev = p;
   }

   return(NULL);
}


