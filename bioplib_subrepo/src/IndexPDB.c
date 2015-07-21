/*************************************************************************

   Program:    
   File:       IndexPDB.c
   
   Version:    V2.0R
   Date:       01.03.94
   Function:   Create an array of pointers into a PDB linked list
   
   Copyright:  (c) SciTech Software 1993-4
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
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
   IndexPDB() creates an array of pointers to each PDB record in a linked
   list. This allows random access to atoms without having to step through
   the PDB linked list.

**************************************************************************

   Usage:
   ======
   pdb.h must be included before using this routine.

   PDB **indx,
       *pdb;
   int natom;

   indx = IndexPDB(pdb, &natom);

**************************************************************************

   Revision History:
   =================
   V1.0  19.07.90 Original
   V1.0a 15.02.91 Corrected comments to match new standard.
   V1.1  01.06.92 ANSIed and documented, FPU condition added
   V2.0  24.02.94 Completely re-written. Note that the calling format
                  has changed!! NOT BACKWARDLY COMPATIBLE!

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdlib.h>

#include "MathType.h"
#include "SysDefs.h"
#include "pdb.h"
#include "macros.h"

/************************************************************************/
/*>PDB **IndexPDB(PDB *pdb, int *natom)
   ------------------------------------
   Input:    PDB   *pdb        Pointer to the start of a PDB linked list.
   Output:   int   *natom      Number of atoms in the PDB linked list.
   Returns:  PDB   **indx      An array of pointers to the PDB records.
                               NULL if unable to allocate memory.

   Creates an array of pointers to PDB from a linked list. This is used
   to allow array style access to items in the linked list:
   e.g. (indx[23])->x will give the x coordinate of the 23rd item

   19.07.90 Original
   01.06.92 ANSIed and documented.
   24.02.94 Re-written. Now allocates and returns the index.
*/
PDB **IndexPDB(PDB *pdb, int *natom)
{
   PDB *p,
       **indx;
   int i=0;

   /* Count the number of entries                                       */
   for(p=pdb, i=0; p!=NULL; NEXT(p)) i++;
   *natom = i;

   /* Allocate memory for the index array                               */
   if((indx = (PDB **)malloc((i+1) * sizeof(PDB *)))==NULL)
      return(NULL);
   

   for(p=pdb, i=0; p!=NULL; NEXT(p))
      indx[i++] = p;

   indx[i] = NULL;
   
   return(indx);
}

