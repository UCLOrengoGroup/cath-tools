/*************************************************************************

   Program:    
   File:       SelAtPDB.c
   
   Version:    V1.8
   Date:       03.20.09
   Function:   Select a subset of atom types from a PDB linked list
   
   Copyright:  (c) SciTech Software 1990-2009
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
   pdbout = selectatoms(pdbin,nsel,sel,natom)

   This routine takes a linked list of type PDB and returns a list
   containing only those atom types specfied in the sel array.

   Input:   pdbin    *PDB      Input list
            nsel     int       Number of atom types to keep
            sel      **char    List of atom types to keep
   Output:  natom    *int      Number of atoms kept
   Returns: pdbout   *PDB      Output list

   To set up the list of atoms to keep, define an array of pointers 
   to char:
   e.g.     char *sel[10]
   Then define the atoms in the list thus:
            SELECT(sel[0],"N   ");
            SELECT(sel[1],"CA  ");
            SELECT(sel[2],"C   ");
            SELECT(sel[3],"O   ");
   The SELECT macro returns a character pointer which will be NULL if
   the allocation it performs fails.

   N.B. The routines are non-destructive; i.e. the original PDB linked 
   list is intact after the selection process

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
   V1.8  04.02.09 SelectAtomsPDB(): Initialize q for fussy compliers

*************************************************************************/
/* Includes
*/
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "SysDefs.h"
#include "MathType.h"

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
/*>PDB *SelectAtomsPDB(PDB *pdbin, int nsel, char **sel, int *natom)
   -----------------------------------------------------------------
   Input:   pdbin    *PDB      Input list
            nsel     int       Number of atom types to keep
            sel      **char    List of atom types to keep
   Output:  natom    *int      Number of atoms kept
   Returns:          *PDB      Output list

   Take a PDB linked list and returns a list containing only those atom 
   types specified in the sel array.

   To set up the list of atoms to keep, define an array of pointers 
   to char:
   e.g.     char *sel[10]
   Then define the atoms in the list thus:
            SELECT(sel[0],"N   ");
            SELECT(sel[1],"CA  ");
            SELECT(sel[2],"C   ");
            SELECT(sel[3],"O   ");
   Ensure the spaces are used!!

   N.B. The routine is non-destructive; i.e. the original PDB linked 
        list is intact after the selection process

   01.03.90 Original    By: ACRM
   28.03.90 Modified to match new version of pdb.h
   24.05.90 Fixed so the variables passed in as sel[] don't 
            *have* to be 4 chars.
   17.05.93 Modified for book. Returns BOOL.
   09.07.93 Modified to return PDB pointer. Changed allocation 
            scheme. Changed back to sel[] variables *must* be 4
            chars.
   04.02.09 Initialize q for fussy compliers
*/
PDB *SelectAtomsPDB(PDB *pdbin, int nsel, char **sel, int *natom)
{
   PDB   *pdbout  = NULL,
         *p,
         *q = NULL;
   int   i;
    
   *natom = 0;
   
   /* Step through the input PDB linked list                            */
   for(p=pdbin; p!= NULL; NEXT(p))
   {
      /* Step through the selection list                                */
      for(i=0; i<nsel; i++)
      {
         /* See if there is a match                                     */
         if(!strncmp(p->atnam,sel[i],4))
         {
            /* Alloacte a new entry                                     */
            if(pdbout==NULL)
            {
               INIT(pdbout, PDB);
               q = pdbout;
            }
            else
            {
               ALLOCNEXT(q, PDB);
            }
            
            /* If failed, free anything allocated and return            */
            if(q==NULL)
            {
               if(pdbout != NULL) FREELIST(pdbout,PDB);
               *natom = 0;
               return(NULL);
            }
            
            /* Increment atom count                                     */
            (*natom)++;
            
            /* Copy the record to the output list (sets ->next to NULL) */
            CopyPDB(q, p);
            
            break;
         }
      }
   }

   /* Return pointer to start of output list                            */
   return(pdbout);
}

