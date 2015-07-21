/*************************************************************************

   Program:    
   File:       OrderPDB.c
   
   Version:    V1.3R
   Date:       23.06.08
   Function:   Functions to modify atom order in PDB linked list
   
   Copyright:  (c) SciTech Software 1993-2008
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
   V1.0  22.02.04 Original
   V1.1  09.03.94 Bug fix in ShuffleResPDB(). Added PCA to sAtoms table
   V1.2  18.03.94 Bug fix in ShuffleResPDB().
   V1.3  22.06.08 Bug fix in ShuffleBB()

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "MathType.h"
#include "macros.h"
#include "pdb.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/
static char sAtoms[MAXSTDAA][MAXATINRES+1][8] = 
{{"ALA","N   ","CA  ","C   ","O   ","CB  ","    ","    ","    ","    ",
        "    ","    ","    ","    ","    "},
 {"ARG","N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","NE  ","CZ  ",
        "NH1 ","NH2 ","    ","    ","    "},
 {"ASN","N   ","CA  ","C   ","O   ","CB  ","CG  ","OD1 ","ND2 ","    ",
        "    ","    ","    ","    ","    "},
 {"ASP","N   ","CA  ","C   ","O   ","CB  ","CG  ","OD1 ","OD2 ","    ",
        "    ","    ","    ","    ","    "},
 {"CYS","N   ","CA  ","C   ","O   ","CB  ","SG  ","    ","    ","    ",
        "    ","    ","    ","    ","    "},
 {"GLN","N   ","CA  ","C   ","O   ","CB  ","CG  ","OE1 ","NE2 ","    ",
        "    ","    ","    ","    ","    "},
 {"GLU","N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","OE1 ","OE2 ",
        "    ","    ","    ","    ","    "},
 {"GLY","N   ","CA  ","C   ","O   ","    ","    ","    ","    ","    ",
        "    ","    ","    ","    ","    "},
 {"HIS","N   ","CA  ","C   ","O   ","CB  ","CG  ","ND1 ","CD2 ","CE1 ",
        "NE2 ","    ","    ","    ","    "},
 {"ILE","N   ","CA  ","C   ","O   ","CB  ","CG1 ","CG2 ","CD  ","    ",
        "    ","    ","    ","    ","    "},
 {"LEU","N   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ","CD2 ","    ",
        "    ","    ","    ","    ","    "},
 {"LYS","N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","CE  ","NZ  ",
        "    ","    ","    ","    ","    "},
 {"MET","N   ","CA  ","C   ","O   ","CB  ","CG  ","SD  ","CE  ","    ",
        "    ","    ","    ","    ","    "},
 {"PHE","N   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ","CD2 ","CE1 ",
        "CE2 ","CZ  ","    ","    ","    "},
 {"PRO","N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","    ","    ",
        "    ","    ","    ","    ","    "},
 {"SER","N   ","CA  ","C   ","O   ","CB  ","OG  ","    ","    ","    ",
        "    ","    ","    ","    ","    "},
 {"THR","N   ","CA  ","C   ","O   ","CB  ","OG1 ","CG2 ","    ","    ",
        "    ","    ","    ","    ","    "},
 {"TRP","N   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ","CD2 ", "NE1 ",
        "CE2 ","CE3 ","CZ2 ","CZ3 ","CH2 "},
 {"TYR","N   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ","CD2 ", "CE1 ",
        "CE2 ","CZ  ","OH  ","    ","    "},
 {"VAL","N   ","CA  ","C   ","O   ","CB  ","CG1 ","CG2 ","    ","    ",
        "    ","    ","    ","    ","    "},
 {"PCA","N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","OE  ","    ",
        "    ","    ","    ","    ","    "}
}  ;


/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>PDB *FixOrderPDB(PDB *pdb, BOOL Pad, BOOL Renum)
   ------------------------------------------------
   Input:   PDB   *pdb    PDB linked list to fix atom order
            BOOL  pad     TRUE: Create dummy coordinate atoms for any
                          missing atoms in standard residues
            BOOL  Renum   TRUE: Renumber the atoms
   Returns: PDB   *       Corrected PDB linked list

   Runs through a PDB linked list and corrects the atom order to match
   the N,CA,C,O,s/c standard. Only standard amino acids are processed.
   The input linked list is modified (i.e. a new list is not built),
   but the return value from the routine should be used for the corrected
   list rather than the input PDB pointer, since the start of the list
   may have changed. i.e. the routine should be called with the form:

                  pdb = FixOrderPDB(pdb,TRUE,TRUE);

   08.07.93 Original    By: ACRM
*/
PDB *FixOrderPDB(PDB *pdb, BOOL Pad, BOOL Renum)
{
   PDB   *start   = NULL,
         *end     = NULL,
         *p       = NULL,
         *ret     = NULL,
         *current = NULL;
   
   for(start=pdb; start!=NULL; start=end)
   {
      /* Find a residue's limits                                        */
      end = FindEndPDB(start);
      
      p = ShuffleResPDB(start, end, Pad);
      
      if(ret == NULL)
      {
         current = ret = p;
      }
      else
      {
         current->next = p;
      }
      
      while(current->next != end)
         NEXT(current);
   }

   if(Renum) RenumAtomsPDB(ret);
   
   return(ret);
}


/************************************************************************/
/*>PDB *ShuffleResPDB(PDB *start, PDB *end, BOOL Pad)
   --------------------------------------------------
   Input:   PDB   *start   Start of residue to be shuffled
            PDB   *end     Start of next residue in linked list (NULL
                           for last residue)
            BOOL  Pad      TRUE: Create dummy records for missing atoms
   Returns: PDB   *        Pointer to new start of linked list

   Shuffle atoms within a residue into the standard order. Returns a 
   pointer to the new first atom in the residue. Atoms not in the known
   list are discarded.
   If we fail to allocate memory for extra atom records, no action is
   taken.

   08.07.93 Original    By: ACRM
   09.03.94 Correctly handles residues not found in the standard list
            (i.e. returns them unmodified)
   17.03.94 If no atoms are found, then we return start. This is
            the case when partial occupancy atoms are named as "N  A",
            "N  B", etc.
*/
PDB *ShuffleResPDB(PDB *start, PDB *end, BOOL Pad)
{
   int   i,
         j;
   char  atnam[8];
   PDB   *ret     = NULL,
         *p       = NULL,
         *extra   = NULL;
   BOOL  found    = FALSE;
         
   /* Search for atom table for this amino acid                         */
   for(i=0; i<MAXSTDAA; i++)
   {
      if(!strncmp(start->resnam,sAtoms[i][0],3))
      {
         /* Found the data for this residue.
            Get each atom name in turn
         */
         for(j=1; j<=MAXATINRES; j++)
         {
            /* Copy the current atom name (for convenience)             */
            strcpy(atnam,sAtoms[i][j]);

            /* Flag for atom found                                      */
            found = FALSE;
            
            /* Break out if we've run out of atom names                 */
            if(atnam[0] == ' ') break;
            
            for(p=start; p!=end; NEXT(p))
            {
               if(!strncmp(p->atnam, atnam, 4) ||
                  (!strncmp(p->resnam,"ILE ",4) &&
                   !strncmp(p->atnam, "CD1 ",4) &&
                   !strncmp(atnam,    "CD  ",4)))
               {
                  /* Atom found, move to return list                    */
                  MovePDB(p, &start, &ret);
                  found = TRUE;
                  break;
               }
            }

            if(Pad && !found)
            {
               INIT(extra,PDB);
               
               if(extra != NULL)
               {
                  /* Copy in information for this residue               */
                  if(ret != NULL)
                     CopyPDB(extra, ret);
                  else
                     CopyPDB(extra, start);
   
                  /* Set required atom name and NULL coordinates        */
                  strcpy(extra->atnam,atnam);
                  extra->x    = (REAL)9999.0;
                  extra->y    = (REAL)9999.0;
                  extra->z    = (REAL)9999.0;
                  extra->occ  = (REAL)0.0;
                  extra->bval = (REAL)20.0;
   
                  /* Now move this record into the return list          */
                  MovePDB(extra, &extra, &ret);
               }
            }
         }
         break;
      } 
   }

   /* If we've searched the whole list without finding our residue type
      simply return the start of the list
   */
   if(i==MAXSTDAA) return(start);

   /* If we found our residue type, but we didn't find any atoms, so
      ret is NULL, again we return the start of the list
   */
   if(ret==NULL) return(start);

   /* Rejoin the shuffled return list to the rest of the list           */
   p=ret;
   LAST(p);
   p->next = end;

   /* Delete any atom records not transferred                           */
   if(start != end)
   {
      /* Move to end of discard list                                    */
      for(p=start; p->next!=end; NEXT(p));
      /* Terminate the list there                                       */
      p->next = NULL;
      /* Free the discard list                                          */
      FREELIST(start, PDB);
   }

   /* Return start of shuffled list                                     */
   return(ret);   
}


/************************************************************************/
/*>BOOL GetAtomTypes(char *resnam, char **AtomTypes)
   -------------------------------------------------
   Input:   char  *resnam      Residue name for which to search
   Output:  char  **AtomTypes  Array of atom names contained in the 
                               residue
   Returns: BOOL               Success

   Fill in atom types for a given residue. AtomTypes must be 
   pre-allocated and must be set up as an array of character pointers
   using Array2D() (or equivalent) rather than being created as a simple
   char types[][] - THIS WILL NOT WORK!

   08.07.93 Original    By: ACRM
*/
BOOL GetAtomTypes(char *resnam, char **AtomTypes)
{
   int   i,
         j;

   for(i=0; i<MAXSTDAA; i++)
   {
      if(!strncmp(resnam,sAtoms[i][0],3))
      {
         for(j=0; j<MAXATINRES; j++)
         {
            strcpy(AtomTypes[j],sAtoms[i][j+1]);
         }
         return(TRUE);
      }
   }

   return(FALSE);
}

/************************************************************************/
/*>PDB *ShuffleBB(PDB *pdb)
   ------------------------
   Input:   PDB  *pdb   Input PDB linked list
   Returns: PDB  *      Start of shuffled PDB linked list

   Shuffles the PDB list to match the standard of N,CA,C,O,CB,other.
   Basically designed to be used with backbones only since the sidechain
   order is not modified. Returns the start of the shuffled PDB linked 
   list.

   13.05.92 Original
   22.06.08 Fixed to use FindNextResidue(). It was continuing to search
            out of the current residue, so if the same residue number in
            a different chain was of the same type, it ended up changing
            that instead!
*/
PDB *ShuffleBB(PDB *pdb)
{
   PDB   *N    = NULL,
         *CA   = NULL,
         *C    = NULL,
         *O    = NULL,
         *CB   = NULL,
         *ret  = NULL,
         *p,
         *next = NULL;

   next = FindNextResidue(pdb);
          
   for(p=pdb; p!=next; NEXT(p))
   {
      if(!strncmp(p->atnam,"N   ",4))
         N  = p;
      else if(!strncmp(p->atnam,"CA  ",4))
         CA = p;
      else if(!strncmp(p->atnam,"C   ",4))
         C  = p;
      else if(!strncmp(p->atnam,"O   ",4))
         O  = p;
      else if(!strncmp(p->atnam,"CB  ",4))
         CB = p;
   }
    
   /* If we didn't find N, just return                                  */
   if(N==NULL) return(pdb);
 
   /* Move atoms in order from the pdb list to the ret list             */
   MovePDB(N,  &pdb, &ret);
   if(CA != NULL) MovePDB(CA, &pdb, &ret);
   if(C  != NULL) MovePDB(C,  &pdb, &ret);
   if(O  != NULL) MovePDB(O,  &pdb, &ret);
   if(CB != NULL) MovePDB(CB, &pdb, &ret);
 
   /* Append the remains of pdb onto ret                                */
   AppendPDB(ret, pdb);
    
   return(ret);
}
