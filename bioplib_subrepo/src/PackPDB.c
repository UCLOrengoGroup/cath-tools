/*************************************************************************

   Program:    
   File:       PackPDB.c
   
   Version:    V1.1R
   Date:       30.05.02
   Function:   Routines for handling packed PDB files
   
   Copyright:  (c) SciTech Software 1993-2002
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  08.07.93 Original
   V1.1  14.03.94 Warning fix in WritePackedResidue()
   V1.2  30.05.02 Changed PDB field from 'junk' to 'record_type'

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "macros.h"
#include "pdb.h"
#include "fsscanf.h"
#include "array.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF      160

/************************************************************************/
/*>BOOL UnPackPDB(FILE *in, FILE *out)
   -----------------------------------
   Input:   FILE  *in    Input file pointer (packed PDB)
            FILE  *out   Output file pointer (normal PDB)
   Returns: BOOL         Success

   Read a packed PDB file and write as a normal (unpacked) file

   08.07.93 Original    By: ACRM
*/
BOOL UnPackPDB(FILE *in, FILE *out)
{
   char  buffer[MAXBUFF];
   int   natom;
   PDB   *pdb     = NULL;
   BOOL  GotHet   = FALSE;
   
   /* Check this is a packed PDB file                                   */
   if(fgets(buffer,MAXBUFF-1,in))
   {
      if(!strncmp(buffer,"PCKPDB",6))
      {
         /* Copy any non-ATOM initial cards                             */
         while(fgets(buffer,MAXBUFF-1,in))
         {
            if(!strncmp(buffer,"HETATM",6) ||
               !strncmp(buffer,"CONECT",6) ||
               !strncmp(buffer,"MASTER",6) ||
               !strncmp(buffer,"END",3))
            {
               GotHet = TRUE;
               break;
            }

            if(!strncmp(buffer,"RES: ", 5))
               break;

            fputs(buffer, out);
         }

         /* Now read the packed ATOM info into a PDB linked list        */
         rewind(in);
         pdb = ReadPackedPDB(in, &natom);

         /* Write the PDB linked list                                   */
         WritePDB(out, pdb);
         
         /* Free the linked list                                        */
         FREELIST(pdb,PDB);
         
         /* Add any records which sit at the end of the file            */
         if(GotHet)
         {
            rewind(in);
            while(fgets(buffer,MAXBUFF-1,in))
            {
               if(!strncmp(buffer,"HETATM",6) ||
                  !strncmp(buffer,"CONECT",6) ||
                  !strncmp(buffer,"MASTER",6) ||
                  !strncmp(buffer,"END",3))
               {
                  fputs(buffer,out);
               }
               
               if(!strncmp(buffer,"RES: ",5))
                  break;
            }
         }
         
         return(TRUE);
      }
   }
   
   return(FALSE);
}
   
   
/************************************************************************/
/*>PDB *ReadPackedPDB(FILE *in, int *natom)
   ----------------------------------------
   Input:   FILE  *in    Input packed PDB file pointer
   Output:  int   *natom Number of atoms read
   Returns: PDB   *      Allocated PDB linked list

   Read a packed PDB file into a PDB linked list. Returns a pointer to
   the start of the list or NULL if failed
   08.07.93 Original    By: ACRM
*/
PDB *ReadPackedPDB(FILE *in, int *natom)
{
   char  buffer[MAXBUFF],
         resnam[8],
         chain[8],
         insert[8],
         **AtomTypes;
   int   resnum   = 0,
         atnum    = 0,
         atcount  = 0;
   REAL  x, y, z, occ, bval;
   PDB   *pdb     = NULL,
         *p       = NULL;
   BOOL  Start    = TRUE;
   
   /* Check this is a packed PDB file                                   */
   if(fgets(buffer,MAXBUFF-1,in))
   {
      if(!strncmp(buffer,"PCKPDB",6))
      {
         /* Allocate a 2D array for the atom types. The first dimension
            must be at least as large as MAXATINRES.
         */
         AtomTypes = (char **)Array2D(sizeof(char),MAXATINRES+1,8);
         
         /* It is a packed PDB file, so unpack into a PDB linked list   */
         while(fgets(buffer,MAXBUFF-1,in))
         {
            if(!strncmp(buffer,"RES:",4))
            {
               TERMINATE(buffer);
               
               /* Read residue information                              */
               fsscanf(buffer,"%5x%4s%1s%4d%1s",
                       resnam,chain,&resnum,insert);
               
               /* Get atom info for this residue type                   */
               GetAtomTypes(resnam, AtomTypes);
               atcount = 0;
               
               /* Clear flag for just copying info                      */
               Start = FALSE;
            }
            else if(!Start)
            {
               TERMINATE(buffer);
               
               /* Read coordinates                                      */
               sscanf(buffer,"%lf %lf %lf %lf %lf",&x,&y,&z,&occ,&bval);
               
               /* Allocate space in linked list                         */
               if(pdb == NULL)
               {
                  INIT(pdb, PDB);
                  p = pdb;
               }
               else
               {
                  ALLOCNEXT(p, PDB);
               }
               
               /* Check space allocated                                 */
               if(p==NULL)
               {
                  if(pdb != NULL) FREELIST(pdb,PDB);
                  return(NULL);
               }
               
               /* Store data in linked list item                        */
               p->x = x;
               p->y = y;
               p->z = z;
               p->occ = occ;
               p->bval = bval;
               p->atnum = ++atnum;
               p->resnum = resnum;
               strcpy(p->record_type, "ATOM  ");
               strcpy(p->atnam,       AtomTypes[atcount++]);
               strcpy(p->resnam,      resnam);
               strcpy(p->insert,      insert);
               strcpy(p->chain,       chain);
               p->next = NULL;
            }
         }
         
         FreeArray2D(AtomTypes,15,8);

         return(pdb);
      }
   }
   
   return(NULL);
}

/************************************************************************/
/*>BOOL PackPDB(FILE *in, FILE *out)
   ---------------------------------
   Input:   FILE  *in    Input normal PDB file pointer
            FILE  *out   Output packed PDB file pointer
   Returns: BOOL         Success

   Pack a PDB file
   08.07.93 Original    By: ACRM
*/
BOOL PackPDB(FILE *in, FILE *out)
{
   PDB   *pdb,
         *p,
         *q;
   int   natom;
   char  buffer[MAXBUFF];

   /* Write the PACK header to the output file                          */
   fprintf(out,"PCKPDB Packed PDB File\n");
   
   /* Copy any non-atom records to the output file                      */
   while(fgets(buffer,MAXBUFF-1,in))
   {
      /* If not ATOM or TER card, just copy it                          */
      if(strncmp(buffer,"ATOM  ",6) && strncmp(buffer,"TER",3))
      {
         fputs(buffer,out);
      }
   }
   
   rewind(in);

   /* Read the PDB file                                                 */
   if((pdb = ReadPDBAtoms(in, &natom)) != NULL)
   {
      pdb = FixOrderPDB(pdb, TRUE, TRUE); /* Pad and renumber           */
      
      for(p=pdb; p!=NULL; p=q)
      {
         q = FindEndPDB(p);
         WritePackedResidue(out, p, q);
      }
   }
   else
   {
      return(FALSE);
   }
   
   if(pdb!=NULL) FREELIST(pdb, PDB);
   
   return(TRUE);
}

/************************************************************************/
/*>void WritePackedResidue(FILE *out, PDB *start, PDB *end)
   --------------------------------------------------------
   Input:   FILE  *out     Output packed PDB file pointer
            PDB   *start   Pointer to start of residue in PDB linked list
            PDB   *end     Pointer to start of next residue in PDB linked
                           list

   Write a PDB residue in packed format
   08.07.93 Original    By: ACRM
   14.03.94 %lf -> %f
*/
void WritePackedResidue(FILE *out, PDB *start, PDB *end)
{
   PDB *p;

   fprintf(out,"RES: %-4s%1s%4d%1s\n", start->resnam, start->chain, 
                                       start->resnum, start->insert);
   for(p=start; p!=end; NEXT(p))
   {
      fprintf(out,"%.3f %.3f %.3f %.2f %.2f\n",
              p->x,
              p->y,
              p->z,
              p->occ,
              p->bval);
   }
}




