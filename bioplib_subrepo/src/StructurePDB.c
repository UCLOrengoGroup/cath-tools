/*************************************************************************

   Program:    
   File:       StructurePDB.c
   
   Version:    V1.1
   Date:       19.05.10
   Function:   Build a structured PDB representation
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2009-10
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

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
   V1.0   24.11.09  Original
   V1.1   19.05.10  PDBRESIDUE and PDBCHAIN are doubly linked lists
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
/*>PDBSTRUCT *AllocPDBStructure(PDB *pdb)
   --------------------------------------
   Input:    PDB        *pdb    PDB linked list
   Returns:  PDBSTRUCT  *       PDBSTRUCT structure containing chain
                                linked list

   Takes a PDB linked list and converts it into a hierarchical structure
   of chains, residues and atoms

   24.11.09  Original   By: ACRM
   19.05.10  Correctly do doubly linked lists
   02.06.10  Fixed some initializations!
*/
PDBSTRUCT *AllocPDBStructure(PDB *pdb)
{
   PDB        *pdbc, *pdbr;
   PDB        *nextchain, *nextres;
   PDBSTRUCT  *pdbstruct = NULL;
   PDBCHAIN   *chain = NULL;
   PDBRESIDUE *residue = NULL;
   
   
   /* Initialization                                                    */
   if((pdbstruct = (PDBSTRUCT *)malloc(sizeof(PDBSTRUCT)))==NULL)
      return(NULL);
   pdbstruct->pdb = pdb;
   pdbstruct->chains = NULL;
   pdbstruct->extras = NULL;             /* 02.06.10                    */
   
   /* Build the chain list                                              */
   for(pdbc=pdb; pdbc!=NULL; pdbc=nextchain)
   {
      nextchain = FindNextChain(pdbc);

      if(pdbstruct->chains == NULL)
      {
         INITPREV(pdbstruct->chains, PDBCHAIN);
         chain = pdbstruct->chains;
      }
      else
      {
         ALLOCNEXTPREV(chain, PDBCHAIN);
      }
      if(chain == NULL)
      {
         FreePDBStructure(pdbstruct);
         return(NULL);
      }
      
      chain->start    = pdbc;
      chain->stop     = nextchain;
      chain->residues = NULL;              /* 02.06.10                  */
      chain->extras   = NULL;              /* 02.06.10                  */
      strcpy(chain->chain, pdbc->chain);
   }
   
   /* Build the residue list for each chain                             */
   for(chain = pdbstruct->chains; chain!=NULL; NEXT(chain))
   {
      for(pdbr = chain->start; pdbr != chain->stop; pdbr=nextres)
      {
         char resid[24];
         int  i, j;
         
         nextres = FindNextResidue(pdbr);
         
         if(chain->residues == NULL)
         {
            INITPREV(chain->residues, PDBRESIDUE);
            residue = chain->residues;
         }
         else
         {
            ALLOCNEXTPREV(residue, PDBRESIDUE);
         }
         if(residue == NULL)
         {
            FreePDBStructure(pdbstruct);
            return(NULL);
         }
         
         residue->start = pdbr;
         residue->stop  = nextres;
         strcpy(residue->chain, pdbr->chain);
         strcpy(residue->insert, pdbr->insert);
         strcpy(residue->resnam, pdbr->resnam);
         residue->resnum = pdbr->resnum;
         residue->extras = NULL;            /* 02.06.10                 */

         /* Build a residue label                                       */
         sprintf(resid, "%s.%d%s", 
                 residue->chain, residue->resnum, residue->insert);
         for(i=0, j=0; resid[i]; i++)
         {
            /* Copy the character if it's not a space, but only copy a 
               full stop if it's not going to be the first character
            */
            if(resid[i] != ' ' && (resid[i] != '.' || j>0))
            {
               residue->resid[j++] = resid[i];
            }
         }
         residue->resid[j++] = '\0';
      }
   }
   
   return(pdbstruct);
}

/************************************************************************/
/*>void FreePDBStructure(PDBSTRUCT *pdbstruct)
   -------------------------------------------
   Input:    PDBSTRUCT  *       PDBSTRUCT structure containing chain
                                linked list

   Frees memory used by the hierarchical description of a PDB structure.
   Note this does not free the underlying PDB linked list

   24.11.09  Original   By: ACRM
*/
void FreePDBStructure(PDBSTRUCT *pdbstruct)
{
   PDBCHAIN *chain = NULL;

   if(pdbstruct == NULL)
      return;
   
   for(chain = pdbstruct->chains; chain!=NULL; NEXT(chain))
   {
      FREELIST(chain->residues, PDBRESIDUE);
   }
   
   FREELIST(pdbstruct->chains, PDBCHAIN);
   free(pdbstruct);
}


/************************************************************************/
/*>PDB *FindNextChain(PDB *pdb)
   ----------------------------
   Input:    PDB        *pdb    PDB linked list
   Returns:  PDB        *       Pointer to start of next chain

   Takes a PDB linked list and find the start of the next chain. This is
   similar to another Bioplib routine which terminates the first chain,
   but this routines doesn't terminate.

   24.11.09  Original   By: ACRM
*/
PDB *FindNextChain(PDB *pdb)
{
   PDB  *p, *ret = NULL;
   char chain;
   
   chain = pdb->chain[0];
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->next == NULL) || (p->next->chain[0] != chain))
      {
         ret = p->next;
         break;
      }
   }

   return(ret);
}

/************************************************************************/
#ifdef DEMO
#include <stdio.h>
int main(int argc, char **argv)
{
   FILE *fp = NULL;
   int natoms;
   PDB *pdb, *p;
   PDBSTRUCT *pdbs;
   PDBCHAIN *chain;
   PDBRESIDUE *residue;
   
   
   fp   = fopen(argv[1], "r");
   pdb  = ReadPDB(fp, &natoms);
   pdbs = AllocPDBStructure(pdb);

   for(chain=pdbs->chains; chain!=NULL; NEXT(chain))
   {
      fprintf(stdout,"COMMENT: Chain %s\n", chain->chain);
      for(residue=chain->residues; residue!=NULL; NEXT(residue))
      {
         fprintf(stdout, "COMMENT: Residue %s\n", residue->resid);
         for(p=residue->start; p!=residue->stop; NEXT(p))
         {
            WritePDBRecord(stdout, p);
         }
      }
   }

   for(chain=pdbs->chains; chain!=NULL; NEXT(chain))
   {
      fprintf(stdout,"COMMENT: Chain %s\n", chain->chain);
      for(p=chain->start; p!=chain->stop; NEXT(p))
      {
         WritePDBRecord(stdout, p);
      }
   }

   FreePDBStructure(pdbs);

   return(0);
}

#endif
   
