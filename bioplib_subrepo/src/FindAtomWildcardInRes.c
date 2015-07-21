/*************************************************************************

   Program:    
   File:       FindAtomWildcardInRes.c
   
   Version:    V1.0
   Date:       27.08.96
   Function:   Find an atom within a residue allowing wild cards
   
   Copyright:  (c) Dr. Andrew C. R. Martin, UCL, 1996
   Author:     Dr. Andrew C. R. Martin
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
   V1.0  27.08.96 Original moved from mutmodel  By: ACRM

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
/*>PDB *FindAtomWildcardInRes(PDB *pdb, char *pattern)
   ---------------------------------------------------
   Input:   PDB   *pdb      Pointer to start of a residue
            char  *pattern  Atom name pattern to find
   Returns: PDB   *         Pointer to requested atom or NULL if not
                            found.

   Finds an atom within the residue given as a PDB pointer. Allows 
   single character wildcards. Thus ?G? maybe used for any atom at the
   gamma position.

   Returns the first atom which matches

   27.08.96 Original   By: ACRM
*/
PDB *FindAtomWildcardInRes(PDB *pdb, char *pattern)
{
   PDB  *p, *pNext;
   int  i;
   BOOL ok;
   
   pNext = FindNextResidue(pdb);
   
   for(p=pdb; p!=pNext; NEXT(p))
   {
      if(!strncmp(p->atnam,pattern,4))
         return(p);
      ok = TRUE;
      for(i=0; i<4; i++)
      {
         if((pattern[i] != p->atnam[i]) && pattern[i] != '?')
         {
            ok = FALSE;
            break;
         }
      }
      if(ok)
         return(p);
   }
   return(NULL);
}



