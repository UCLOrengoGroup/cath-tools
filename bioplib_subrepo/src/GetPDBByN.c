/*************************************************************************

   Program:    
   File:       GetPDBByN.c
   
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
/*>PDB *GetPDBByN(PDB *pdb, int n)
   -------------------------------
   Input:   PDB   *pdb    PDB linked list
            int   n       Offset into linked list
   Returns: PDB   *       Pointer to n'th item in linked list

   Gets a pointer to a pdb item by taking a PDB linked list and an 
   integer.
   The pointer returned is the n'th item in the list
   
   13.05.92 Original
*/
PDB *GetPDBByN(PDB *pdb,
               int n)
{
   PDB *p;
   int i;
   
   for(i=0, p=pdb; p && i<n ; NEXT(p), i++) ;

   return(p);
}

