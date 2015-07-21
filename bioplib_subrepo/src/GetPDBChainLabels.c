/*************************************************************************

   Program:    
   File:       GetPDBChainLabels.c
   
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
/*>char *GetPDBChainLabels(PDB *pdb)
   ---------------------------------
   Input:   PDB    *pdb      PDB linked list
   Returns: char   *         Allocated string containing chain labels
                             NULL if unable to allocate memory

   Scans a PDB linked list for chain names. Allocates memory for a 
   string containing these labels which is returned.

   N.B. You must free the allocated memory when you've finished with it!

   25.07.95 Original    By: ACRM
*/
char *GetPDBChainLabels(PDB *pdb)
{
   char *chains;
   int  nchains   = 0,
        maxchains = 16;
   PDB  *p;
   
   /* Just return if linked list is NULL                                */
   if(pdb==NULL)
      return(NULL);

   /* Allocate a chunk for storing the chains                           */
   if((chains = (char *)malloc(maxchains * sizeof(char)))==NULL)
      return(NULL);

   /* Set up first chain label                                          */
   chains[nchains] = pdb->chain[0];

   /* Run through the linked list                                       */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If chain label has changed                                     */
      if(p->chain[0] != chains[nchains])
      {
         /* Increment chain count and reallocate memory if needed       */
         if(++nchains == maxchains)
         {
            maxchains += 16;
            if((chains = realloc(chains, maxchains * sizeof(char)))==NULL)
               return(NULL);
         }
         /* Store this new chain label                                  */
         chains[nchains] = p->chain[0];
      }
   }

   /* Increment chain count and reallocate memory if needed             */
   if(++nchains == maxchains)
   {
      maxchains += 16;
      if((chains = realloc(chains, maxchains * sizeof(char)))==NULL)
         return(NULL);
   }

   /* Terminate the chain list with a NUL character                     */
   chains[nchains] = '\0';

   return(chains);
}


