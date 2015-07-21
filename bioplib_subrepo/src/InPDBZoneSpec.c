/*************************************************************************

   Program:    
   File:       InPDBZoneSpec.c
   
   Version:    V1.3
   Date:       19.09.96
   Function:   
   
   Copyright:  (c) SciTech Software 1993-6
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
   V1.0  30.09.92 Original
   V1.1  16.06.93 Tidied for book. Mode now a char.
   V1.2  18.06.96 Added InPDBZone() from QTree program
   V1.3  19.09.96 Added InPDBZoneSpec()

*************************************************************************/
/* Includes
*/
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
/*>BOOL InPDBZoneSpec(PDB *p, char *resspec1, char *resspec2)
   ----------------------------------------------------------
   Input:   PDB    *p         Pointer to a PDB record
            char   *resspec1  Res spec for first residue
            char   *resspec2  Res spec for last residue
   Returns: BOOL              Is p in the range specified?

   Determines whether a PDB pointer is within a residue range specified
   using standard format: [c]nnn[i]

   Also handles the residue spec of c* (i.e. chain name and a * to
   indicate all residues in a given chain). This must be given as
   resspec1 (resspec2 is then ignored).

   Calls InPDBZone() to do the actual work

   19.09.96 Original  By: ACRM
*/
BOOL InPDBZoneSpec(PDB *p, char *resspec1, char *resspec2)
{
   char chain1,  chain2,
        insert1, insert2;
   int  res1,    res2;

   /* Check for wildcard specification of whole chain                   */
   if(resspec1[1] == '*')
   {
      UPPER(resspec1);
      if(p->chain[0] == resspec1[0])
      {
         return(TRUE);
      }
      else
      {
         return(FALSE);
      }
   }
   
   ParseResSpec(resspec1, &chain1, &res1, &insert1);
   ParseResSpec(resspec2, &chain2, &res2, &insert2);

   if(chain1 != chain2)
      return(FALSE);
   
   return(InPDBZone(p, chain1, res1, insert1, res2, insert2));
}
