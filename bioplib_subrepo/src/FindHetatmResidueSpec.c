/*************************************************************************

   Program:    
   File:       FindHetatmResidueSpec.c
   
   Version:    V1.0
   Date:       26.10.11
   Function:   Parse a residue specification
   
   Copyright:  (c) 
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
   V1.0  26.10.11 Original based on FindResidueSpec.c

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
/*>PDB *FindHetatmResidueSpec(PDB *pdb, char *resspec)
   ---------------------------------------------------
   Input:   PDB   *pdb      PDB linked list
            char  *resspec  Residue specification
   Returns: PDB   *         Pointer to first atom of specified residue
                            (NULL if not found).

   Search a PDB linked list for a specified residue (given as
   [chain]num[insert]) but limits search to HETATM residues

   26.10.11 Original    By: ACRM
*/
PDB *FindHetatmResidueSpec(PDB *pdb, char *resspec)
{
   char chain,
        insert;
   int  resnum;

   if(ParseResSpec(resspec, &chain, &resnum, &insert))    
      return(FindHetatmResidue(pdb, chain, resnum, insert));
   
   return(NULL);
}

