/*************************************************************************

   Program:    
   File:       RenumAtomsPDB.c
   
   Version:    V1.2
   Date:       27.02.98
   Function:   
   
   Copyright:  (c) SciTech Software 1993-8
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
   V1.1  01.03.94 Original
   V1.2  27.02.98 Removed unreachable break from switch()

*************************************************************************/
/* Includes
*/
#include "macros.h"
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
/*>void RenumAtomsPDB(PDB *pdb)
   ----------------------------
   I/O:   PDB  *pdb   PDB linked list to renumber

   Renumber the atoms throughout a PDB linked list
   01.08.93 Original
*/
void RenumAtomsPDB(PDB *pdb)
{
   PDB *p;
   int i;

   for(p=pdb, i=1; p!=NULL; NEXT(p)) 
      p->atnum=i++;
}

