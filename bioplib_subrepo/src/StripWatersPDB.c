/*************************************************************************

   Program:    
   File:       StripWatersPDB.c
   
   Version:    V1.0
   Date:       30.04.08
   Function:   
   
   Copyright:  (c) Dr. Andrew C. R. Martin, UCL, 2008
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
   V1.0  30.04.08 Original based on StripHPDB()   By: ACRM

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
/*>PDB *StripWatersPDB(PDB *pdbin, int *natom)
   -------------------------------------------
   Input:   pdbin    *PDB      Input list
   Output:  natom    *int      Number of atoms kept
   Returns:          *PDB      Output list

   Take a PDB linked list and returns the PDB list minus waters

   N.B. The routine is non-destructive; i.e. the original PDB linked 
        list is intact after the selection process

   30.04.08 Original based on StripHPDB()   By: ACRM
*/
PDB *StripWatersPDB(PDB *pdbin, int *natom)
{
   PDB   *pdbout  = NULL,
         *p,
         *q;
    
   *natom = 0;
   
   /* Step through the input PDB linked list                            */
   for(p=pdbin; p!= NULL; NEXT(p))
   {
      if(!ISWATER(p))
      {
         /* Allocate a new entry                                        */
         if(pdbout==NULL)
         {
            INIT(pdbout, PDB);
            q = pdbout;
         }
         else
         {
            ALLOCNEXT(q, PDB);
         }
         
         /* If failed, free anything allocated and return               */
         if(q==NULL)
         {
            if(pdbout != NULL) FREELIST(pdbout,PDB);
            *natom = 0;
            return(NULL);
         }
         
         /* Increment atom count                                        */
         (*natom)++;
         
         /* Copy the record to the output list (sets ->next to NULL)    */
         CopyPDB(q, p);
      }
   }

   /* Return pointer to start of output list                            */
   return(pdbout);
}
