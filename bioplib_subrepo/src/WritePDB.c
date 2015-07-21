/*************************************************************************

   Program:    
   File:       WritePDB.c
   
   Version:    V1.9R
   Date:       22.09.05
   Function:   Write a PDB file from a linked list
   
   Copyright:  (c) SciTech Software 1993-2005
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
   This routine will write a .PDB file of any size from a linked list of 
   the protein structure. This list is contained in a linked set of 
   structures of type pdb_entry. The strucure is set up by including the 
   file "pdb.h". For details of the structure, see this file.

**************************************************************************

   Usage:
   ======
   WritePDB(fp, pdb)
   Input:   FILE  *fp      A pointer to the file to write
            PDB   *pdb     The start of the PDB linked list.

**************************************************************************

   Revision History:
   =================
   V1.0  08.03.89 Original
   V1.2  28.03.90 Modified to match the correct column definition of 
                  ReadPDB V1.2   (N.B. There was no V1.1)
   V1.3  01.06.92 Corrected header, to match standard. Autodoc'd, 
                  ANSIed. Added FPU check.
   V1.4  10.06.93 Changed to use NEXT() macro. void types
   V1.5  22.02.94 Added TER card at end of file
   V1.6  15.02.01 Writes using atnam_raw so atom name is unchanged from
                  input
   V1.7  30.05.02 Changed PDB field from 'junk' to 'record_type'
   V1.8  03.06.05 'atnam_raw' no longer includes the alternate indicator
                  which is now in 'altpos'
   V1.9  22.09.06 Added WritePDBRecordAtnam()

*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "MathType.h"
#include "pdb.h"
#include "macros.h"

/************************************************************************/
/*>void WritePDB(FILE *fp, PDB *pdb)
   ---------------------------------
   Input:   FILE *fp   PDB file pointer to be written
            PDB  *pdb  PDB linked list to write

   Write a PDB linked list by calls to WritePDBRecord()

   08.03.89 Original
   01.06.92 ANSIed and autodoc'd
   10.06.93 Uses NEXT macro; void type
   08.07.93 Added insertion of TER cards
   22.02.94 And a TER card at the end of the file
*/
void WritePDB(FILE *fp,
              PDB  *pdb)
{
   PDB   *p;
   char  PrevChain[8];
   
   strcpy(PrevChain,pdb->chain);

   for(p = pdb ; p ; NEXT(p))
   {
      if(strncmp(PrevChain,p->chain,1))
      {
         /* Chain change, insert TER card                               */
         fprintf(fp,"TER   \n");
         strcpy(PrevChain,p->chain);
      }
      WritePDBRecord(fp,p);
   }
   fprintf(fp,"TER   \n");
}


/************************************************************************/
/*>void WritePDBRecord(FILE *fp, PDB *pdb)
   ---------------------------------------
   Input:   FILE  *fp     PDB file pointer to be written
            PDB   *pdb    PDB linked list record to write

   Write a PDB record

   08.03.89 Original
   28.03.90 Changed to match ReadPDB() V1.2 for column widths
   01.06.92 ANSIed and autodoc'd
   10.06.93 void type
   22.06.93 Changed to %lf. Ljust strings
   11.03.94 %lf back to %f (!)
   15.02.01 Modified to use atnam_raw
   03.06.05 Modified to use altpos
*/
void WritePDBRecord(FILE *fp,
                    PDB  *pdb)
{
   fprintf(fp,"%-6s%5d %-4s%c%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           pdb->record_type,
           pdb->atnum,
           pdb->atnam_raw,
           pdb->altpos,
           pdb->resnam,
           pdb->chain,
           pdb->resnum,
           pdb->insert,
           pdb->x,
           pdb->y,
           pdb->z,
           pdb->occ,
           pdb->bval);
}
/************************************************************************/
/*>void WritePDBRecordAtnam(FILE *fp, PDB *pdb)
   --------------------------------------------
   Input:   FILE  *fp     PDB file pointer to be written
            PDB   *pdb    PDB linked list record to write

   Write a PDB record

   08.03.89 Original
   28.03.90 Changed to match ReadPDB() V1.2 for column widths
   01.06.92 ANSIed and autodoc'd
   10.06.93 void type
   22.06.93 Changed to %lf. Ljust strings
   11.03.94 %lf back to %f (!)
   15.02.01 Modified to use atnam_raw
   03.06.05 Modified to use altpos
   22.09.05 This is like the old version which used atnam rather
            than atnam_raw
*/
void WritePDBRecordAtnam(FILE *fp,
                         PDB  *pdb)
{
   fprintf(fp,"%-6s%5d  %-4s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           pdb->record_type,
           pdb->atnum,
           pdb->atnam,
           pdb->resnam,
           pdb->chain,
           pdb->resnum,
           pdb->insert,
           pdb->x,
           pdb->y,
           pdb->z,
           pdb->occ,
           pdb->bval);
}
