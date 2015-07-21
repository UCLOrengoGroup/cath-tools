/*************************************************************************

   Program:    
   File:       FindHetatmResidue.c
   
   Version:    V1.0
   Date:       26.10.11
   Function:   Parse a residue specification
   
   Copyright:  
   Author:     Dr. Andrew C. R. Martin
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
   V1.0  26.10.11 Original based on FindResidue.c

*************************************************************************/
/* Includes
*/
#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "pdb.h"
#include "SysDefs.h"
#include "MathType.h"
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
/*>PDB *FindHetatmResidue(PDB *pdb, char chain, int resnum, char insert)
  ----------------------------------------------------------------------
  Finds a pointer to the start of a residue in a PDB linked list, but
  requires the residue is a HETATM record

  26.10.11 Original   By: ACRM
*/
PDB *FindHetatmResidue(PDB *pdb, char chain, int resnum, char insert)
{
   PDB *p;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((!strncmp(p->record_type,"HETATM",6)) &&
      	 (p->resnum    == resnum) &&
         (p->insert[0] == insert) &&
         (p->chain[0]  == chain))
	 {
           return(p);
	 }
   }
   return(NULL);
}


