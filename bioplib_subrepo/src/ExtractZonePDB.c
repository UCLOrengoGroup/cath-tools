/*************************************************************************

   Program:    
   File:       ExtractZonePDB.c
   
   Version:    V1.13R
   Date:       29.10.10
   Function:   PDB linked list manipulation
   
   Copyright:  (c) SciTech Software 1992-2010
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
   V1.11 22.03.05 Extracted range is limited by specified residues
   V1.12 22.03.06 Modified ExtractZonePDB() to allow non-exact ranges
   V1.13 29.10.10 Fixed bug when end of zone was last residue in a chain

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdlib.h>

#include "MathType.h"
#include "SysDefs.h"
#include "pdb.h"
#include "macros.h"
#include "general.h"

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
/*>PDB *ExtractZonePDB(PDB *inpdb, char *chain1, int resnum1, 
                       char *insert1, char *chain2, int resnum2, 
                       char *insert2)
   -----------------------------------------------------------------------
   Input:   PDB    *inpdb   Input PDB linked list
            char   *chain1  Start residue chain name
            int    resnum1  Start residue number
            char   *insert1 Start residue insert code
            char   *chain2  End residue chain name
            int    resnum2  End residue number
            char   *insert2 End residue insert code
   Returns: PDB    *        PDB linked list of the region of interest.

   Reduces a PDB linked list to those residues within a specified zone.
   Note that the PDB linked list is duplicated before extraction so
   pointers do not match those in the input PDB linked list. Excess
   records in the new PDB linked list are freed.

   10.01.96 Original   By: ACRM
   22.03.06 Modified to allow non-exact zones. i.e. the extracted zone
            will be the widest subset of the specified zone. So, if
            you specifiy 30-35Z and the PDB file only has 30-35B
            then that will be extracted.
   29.10.10 Fixed extraction where end of zone matched last residue in
            a chain
*/
PDB *ExtractZonePDB(PDB *inpdb, char *chain1, int resnum1, char *insert1,
                    char *chain2, int resnum2, char *insert2)
{
   PDB *pdb, *p, 
       *start = NULL, 
       *last  = NULL, 
       *prev  = NULL;

   /* Duplicate the PDB linked list                                     */
   if((pdb = DupePDB(inpdb))==NULL)
      return(NULL);

   /* Find the first residue in the PDB linked list                     
      prev will point to the last atom before the first atom in the zone
      start will point to the first atom in the zone
    */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->chain[0] == chain1[0]) &&
         ((p->resnum > resnum1) ||
          ((p->resnum == resnum1) &&
           (p->insert[0] >= insert1[0]))))
      {
         start = p;
         break;
      }
      prev = p;
   }

   if(start==NULL)
   {
      FREELIST(pdb, PDB);
      return(NULL);
   }

   /* Find the last residue in the PDB linked list                      
      last will point to the last atom in the zone

      29.10.10 Also breaks out if chain1 and chain2 are the same but
      we've now come to a different chain. This fixes a bug where
      the code wouldn't break out if resnum2 was the last residue in
      a chain By: ACRM
    */
   for(p=start; p!=NULL; NEXT(p))
   {
      if((p->chain[0] == chain2[0]) && /* If chain is the same and...  */
         ((p->resnum > resnum2) ||     /* Residue number exceeded or.. */
          ((p->resnum == resnum2) &&   /* Resnum same, insert exceeded */
           (p->insert[0] > insert2[0]))))
      {
         break;
      }
      /* Both zone ends are in the same chain so, if we got here we have 
         found the right chain. If the current chain is now different 
         from chain2, then we've gone off the end of the chain
      */
      if((chain1[0]   == chain2[0]) &&
         (p->chain[0] != chain2[0]))
      {
         break;
      }
      last = p;
   }

   if(last==NULL)
   {
      FREELIST(pdb, PDB);
      return(NULL);
   }

   /* Free linked list after 'last'                                     */
   if(last->next != NULL)
   {
      FREELIST(last->next, PDB);
      last->next = NULL;
   }
   
   /* Unlink 'start' from rest of linked list and free memory before 
      'start'
   */
   if(prev != NULL)
   {
      prev->next = NULL;
      FREELIST(pdb, PDB);
   }

   return(start);
}


