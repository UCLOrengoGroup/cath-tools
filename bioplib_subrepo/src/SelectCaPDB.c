/*************************************************************************

   Program:    
   File:       SelectCaPDB.c
   
   Version:    V1.7
   Date:       16.10.96
   Function:   
   
   Copyright:  (c) SciTech Software 1990-6
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
   V1.0  01.03.90 Original   By: ACRM
   V1.1  28.03.90 Modified to match new version of pdb.h
   V1.2  24.05.90 Fixed so the variables passed in as sel[] don't 
                  *have* to be 4 chars.
   V1.3  17.05.93 Modified for book. Returns BOOL.
   V1.4  09.07.93 Modified to return PDB pointer. Changed allocation 
                  scheme. Changed back to sel[] variables *must* be 4
                  chars.
   V1.5  01.11.94 Added HStripPDB()
   V1.6  26.07.95 Removed unused variables
   V1.7  16.10.96 Added SelectCaPDB()

*************************************************************************/
/* Includes
*/
#include <string.h>
#include <stdlib.h>

#include "SysDefs.h"
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
/*>PDB *SelectCaPDB(PDB *pdb)
   --------------------------
   Input:   PDB  *pdb     A PDB linked list
   Output:  PDB  *        The linked list reduced to CA atoms only

   Reduce a PDB linked list to CA atoms only. Returns the reduced PDB
   linked list. Note that the input linked list is modified (storage for
   non-CA atoms is freed) and that the pointer input to the routine
   is likely to become invalid since the first atom is generally no a 
   CA and the first record will therefore be freed.

   Normally you would call this routine with something like:

   pdb = SelectCaPDB(pdb);

   i.e. the same variable name is used for the input and the return.

   This is used in place of the normal SelectAtomsPDB() since it is
   simpler and faster to call for this rather common requirement.

   N.B. Unlike SelectAtomsPDB(), the routine is destructive; i.e. 
        the original PDB linked list is NOT intact after the 
        selection process.

   15.10.96 Original   By: ACRM
   11.01.02 Added check on pdb!=NULL in case there are no CAs (DNA etc)
*/
PDB *SelectCaPDB(PDB *pdb)
{
   PDB *p, *prev;

   /* Remove atoms up to the first CA                                   */
   while((pdb!=NULL) && strncmp(pdb->atnam,"CA  ",4))
   {
      p=pdb->next;
      free(pdb);
      pdb=p;
   }

   /* Remove the rest of the non-CA atoms                               */
   for(p=pdb, prev=pdb; p!=NULL; NEXT(p))
   {
      if(strncmp(p->atnam,"CA  ",4))
      {
         prev->next = p->next;
         free(p);
         p = prev;
      }
      prev = p;
   }

   return(pdb);
}

