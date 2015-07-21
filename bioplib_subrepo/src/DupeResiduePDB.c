/*************************************************************************

   Program:    
   File:       DupeResiduePDB.c
   
   Version:    V1.1
   Date:       08.11.07
   Function:   Create a new PDB linked list with a copy of a residue
   
   Copyright:  (c) Dr. Andrew C. R. Martin, UCL, 1996-2007
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
   V1.0 27.08.96 Original from mutmodel  By: ACRM
   V1.1 08.11.07 Initialize p and q; Moved into bioplib

*************************************************************************/
/* Includes
*/
#include <stdlib.h>
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
/*>PDB *DupeResiduePDB(PDB *in)
   ----------------------------
   Input:   PDB  *in     PDB linked list pointing to residue to duplicate
   Returns: PDB  *       Duplicate linked list of residue at `in'
                         (NULL if allocation fails)

   Makes a duplicate PDB linked list of just the residue pointed to by
   `in'

   27.08.96 Original   By: ACRM
   08.11.07 Initialize p and q
            Moved into bioplib
*/
PDB *DupeResiduePDB(PDB *in)
{
   PDB *pNext,
       *out = NULL,
       *p = NULL, 
       *q = NULL;
   
   /* Find the next residue                                             */
   pNext = FindNextResidue(in);

   for(p=in; p!=pNext; NEXT(p))
   {
      if(out==NULL)
      {
         INIT(out, PDB);
         q=out;
      }
      else
      {
         ALLOCNEXT(q, PDB);
      }
      if(q==NULL)
      {
         FREELIST(out, PDB);
         return(NULL);
      }
      
      CopyPDB(q, p);
   }
   
   return(out);
}


