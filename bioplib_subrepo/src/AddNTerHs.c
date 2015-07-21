/*************************************************************************

   Program:    
   File:       AddNTerHs.c
   
   Version:    V1.6R
   Date:       03.06.05
   Function:   Routines to add N-terminal hydrogens and C-terminal
               oxygens.
   
   Copyright:  (c) SciTech Software 1994-2006
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

   int AddNTerHs(PDB **ppdb, BOOL Charmm)
   --------------------------------------
   This routine is used to generate a set of N-terminal hydrogens.
   By default GROMOS naming is used; the Charmm flag changes this to
   Charmm format in which case an extra NTER residue will be created.
   Any pre-existing backbone H will be removed. Typically this
   routine would be called after calling HAddPDB()

**************************************************************************

   Revision History:
   =================
   V1.0  24.08.94 Original    By: ACRM
   V1.1  05.10.94 Removed unused variables
   V1.2  12.11.96 If any of the antecedant coordinates are undefined, set
                  the terminal oxygen to NULL coordinates
   V1.3  13.11.96 Also checks for missing CA,C and O1 records
   V1.4  30.05.02 Changed PDB field from 'junk' to 'record_type'
   V1.5  06.02.03 Handles atnam_raw
   V1.6  03.06.05 Handles altpos

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SysDefs.h"
#include "MathType.h"
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
static PDB *KillNTerH(PDB *pdb);
static int doAddGromosNTer(PDB **ppdb, PDB *nter);
static int doAddCharmmNTer(PDB **ppdb, PDB *nter);

/************************************************************************/
/*>int AddNTerHs(PDB **ppdb, BOOL Charmm)
   --------------------------------------
   I/O:     PDB   **pdb      Pointer to pointer to PDB linked list
   Input:   BOOL  Charmm     Do Charmm style Nter
   Returns: int              Number of hydrogens added

   Adds hydrogens onto the N-termini

   23.08.94 Original    By: ACRM
*/
int AddNTerHs(PDB **ppdb, BOOL Charmm)
{
   PDB  *p;
   char chain = '-';
   int  nhyd  = 0;
   
   for(p = *ppdb; p!=NULL; NEXT(p))
   {
      if(p->chain[0] != chain)
      {
         chain = p->chain[0];

         /* Kill the H if there is one. Makes the (safe...) assumption 
            that the H won't be the first atom (i.e. ignore the return)
         */
         KillNTerH(p);

         if(Charmm)
            nhyd += doAddCharmmNTer(ppdb, p);
         else
            nhyd += doAddGromosNTer(ppdb, p);
      }
   }
   
   return(nhyd);
}

/************************************************************************/
/*>static PDB *KillNTerH(PDB *pdb)
   -------------------------------
   I/O:     PDB    *pdb      PDB linked list
   Returns: PDB    *         New start of linked list

   Remove the backbone hydrogen from Nter residue

   23.08.94 Original    By: ACRM
*/
static PDB *KillNTerH(PDB *pdb)
{
   PDB *p,
       *end,
       *prev = NULL;

   end = FindEndPDB(pdb);
   
   for(p=pdb; p!=end; NEXT(p))
   {
      if(!strncmp(p->atnam,"H   ",4))
      {
         if(prev==NULL)         /* Is the first item in the list        */
         {
            prev = p->next;
            free(p);
            return(prev);
         }
         else                   /* Is NOT the first item in the list    */
         {
            prev->next = p->next;
            free(p);
            return(pdb);
         }
      }
      prev = p;
   }
   return(pdb);
}
/************************************************************************/
/*>static int doAddGromosNTer(PDB **ppdb, PDB *nter)
   -------------------------------------------------
   I/O:     PDB   **pdb      Pointer to pointer to PDB linked list
   Input:   PDB   *nter      Pointer to an N-terminus
   Returns: int              Number of Hs added (0 if error)

   Does the actual work of adding the hydrogens onto an N terminus

   23.08.94 Original    By: ACRM
   24.08.94 Set B-values of Hs to 20.0
   05.10.94 Removed unused variables
   06.02.03 Handles atnam_raw
   03.06.05 Handles altpos
*/
static int doAddGromosNTer(PDB **ppdb, PDB *nter)
{
   COOR coor[3];
   PDB  *prev = NULL,
        *H1, *H2, *H3;

   /* If not start of list find pointer to previous record              */
   if(*ppdb != nter)
   {
      for(prev = *ppdb; prev->next != nter; NEXT(prev));
   }

   if(CalcTetraHCoords(nter, coor))
   {
      /* Initialise some space to store the 3 extra hydrogens           */
      INIT(H1, PDB);
      INIT(H2, PDB);
      INIT(H3, PDB);
      if(H1==NULL || H2==NULL || H3==NULL) return(0);
      
      /* Initialise the hydrogens with the res info, etc.               */
      CopyPDB(H1, nter);
      CopyPDB(H2, nter);
      CopyPDB(H3, nter);
      H1->bval = H2->bval = H3->bval = (REAL)20.0;

      /* Link these extra records into the linked list                  */
      H3->next   = nter->next;
      H1->next   = H2;
      H2->next   = nter;
      nter->next = H3;

      /* Set the coordinates                                            */
      H1->x = coor[0].x;
      H1->y = coor[0].y;
      H1->z = coor[0].z;
      strcpy(H1->atnam, "H1  ");
      strcpy(H1->atnam_raw, " H1 ");
      H1->altpos = ' ';                   /* 03.06.05                   */
      
      H2->x = coor[1].x;
      H2->y = coor[1].y;
      H2->z = coor[1].z;
      strcpy(H2->atnam, "H2  ");
      strcpy(H2->atnam_raw, " H2 ");
      H2->altpos = ' ';                   /* 03.06.05                   */

      H3->x = coor[2].x;
      H3->y = coor[2].y;
      H3->z = coor[2].z;
      strcpy(H3->atnam, "H3  ");
      strcpy(H3->atnam_raw, " H3 ");
      H3->altpos = ' ';                   /* 03.06.05                   */

      /* Correctly link the new start into the whole linked list        */
      if(prev != NULL)
         prev->next = H1;
      else
         *ppdb = H1;
         
      return(3);
   }
   return(0);
}

/************************************************************************/
/*>static int doAddCharmmNTer(PDB **ppdb, PDB *nter)
   -------------------------------------------------
   I/O:     PDB   **pdb      Pointer to pointer to PDB linked list
   Input:   PDB   *nter      Pointer to an N-terminus
   Returns: int              Number of Hs added (0 if error)

   Does the actual work of adding the hydrogens onto an N terminus

   23.08.94 Original    By: ACRM
   24.08.94 Set B-values of H3 to 20.0
   05.10.94 Removed unused variables
   06.02.03 Handles atnam_raw
   03.06.05 Handles altpos
*/
static int doAddCharmmNTer(PDB **ppdb, PDB *nter)
{
   COOR coor[3];
   PDB  *prev = NULL,
        *H1, *H2, *H3;

   /* If not start of list find pointer to previous record              */
   if(*ppdb != nter)
   {
      for(prev = *ppdb; prev->next != nter; NEXT(prev));
   }

   if(CalcTetraHCoords(nter, coor))
   {
      /* Initialise some space to store the 3 extra hydrogens           */
      INIT(H1, PDB);
      INIT(H2, PDB);
      INIT(H3, PDB);
      if(H1==NULL || H2==NULL || H3==NULL) return(0);
      
      /* Initialise the hydrogens with the res info, etc.               */
      H1->atnum  = 9998;
      H1->resnum = 9999;
      H1->occ    = 1.0;
      H1->bval   = 20.0;
      strcpy(H1->resnam,      "NTER");
      strcpy(H1->atnam,       "HT1 ");
      strcpy(H1->atnam_raw,   " HT1");
      strcpy(H1->record_type, "ATOM  ");
      strcpy(H1->chain,       nter->chain);
      strcpy(H1->insert,      " ");
      H1->altpos = ' ';                   /* 03.06.05                   */
      
      H2->atnum  = 9999;
      H2->resnum = 9999;
      H2->occ    = 1.0;
      H2->bval   = 20.0;
      strcpy(H2->resnam,      "NTER");
      strcpy(H2->atnam,       "HT2 ");
      strcpy(H2->atnam_raw,   " HT2");
      strcpy(H2->record_type, "ATOM  ");
      strcpy(H2->chain,       nter->chain);
      strcpy(H2->insert,      " ");
      H2->altpos = ' ';                   /* 03.06.05                   */

      CopyPDB(H3, nter);
      strcpy(H3->atnam,  "HT3 ");
      strcpy(H3->atnam_raw,   " HT3");
      H3->bval   = 20.0;
      H3->altpos = ' ';                   /* 03.06.05                   */

      /* Link these extra records into the linked list                  */
      H3->next   = nter->next;
      H1->next   = H2;
      H2->next   = nter;
      nter->next = H3;

      /* Set the coordinates                                            */
      H1->x = coor[0].x;
      H1->y = coor[0].y;
      H1->z = coor[0].z;

      H2->x = coor[1].x;
      H2->y = coor[1].y;
      H2->z = coor[1].z;

      H3->x = coor[2].x;
      H3->y = coor[2].y;
      H3->z = coor[2].z;

      /* Change the name of the Nter nitrogen                           */
      strcpy(nter->atnam, "NT  ");
      strcpy(nter->atnam_raw, " NT ");

      /* Correctly link the new start into the whole linked list        */
      if(prev != NULL)
         prev->next = H1;
      else
         *ppdb = H1;
         
      return(3);
   }
   return(0);
}
