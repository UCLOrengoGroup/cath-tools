/*************************************************************************

   Program:    
   File:       ReadPDB.c
   
   Version:    V2.22
   Date:       21.12.11
   Function:   Read coordinates from a PDB file 
   
   Copyright:  (c) SciTech Software 1988-2011
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
   pdb = ReadPDB(fp,natom) - This subroutine will read a .PDB file
   of any size and form a linked list of the protein structure.
   This list is contained in a linked set of structures of type
   pdb_entry. The strucure is set up by including the file
   "pdb.h". For details of the structure, see this file.

   To free the space created by this routine, call FREELIST(pdb,PDB).

   The parameters passed to the subroutine are:
   fp    - A pointer to type FILE in which the .PDB file is stored.
   pdb   - A pointer to type PDB.
   natom - A pointer to type integer in which the number of atoms
           found is stored.

   As of V2.3, the routine makes provision for partial occupancies. If 
   the occupancies are 1.0 or 0.0, the atoms are read verbatim. If not,
   only the highest occupancy atoms are read and the atom names are 
   corrected to remove alternative labels. This behaviour can be 
   overridden by calling one of the ...OccRank() routines to read lower 
   occupancy atoms. If any partial occupancy atoms are read the global
   flag gPDBPartialOcc is set to TRUE.

NOTE:  Although some of the fields are represented by a single character,
       they are still stored in character arrays.

BUGS:  The subroutine cannot read files with VAX Fortran carriage control!
       It just sits there and page faults like crazy.

BUGS:  The multiple occupancy code assumes that all positions for a given
       atom in consecutive records of the file

BUGS:  25.01.05 Note the multiple occupancy code won't work properly for
       3pga where atoms have occupancies of zero and one

**************************************************************************

   Usage:
   ======
   pdb = ReadPDB(fp,natom)
   
   Input:   FILE     *fp      A pointer to type FILE in which the
                              .PDB file is stored.
   Returns: PDB      *pdb     A pointer to the first allocated item of
                              the PDB linked list
   Output:  int      *natom   Number of atoms read.

**************************************************************************

   Revision History:
   =================
   V1.0  04.11.88 Original
   V1.1  07.02.89 Now ignores any records from the .PDB file which 
                  don't start with ATOM or HETATM.
   V1.2  28.03.90 Some fields altered to match the exact specifications 
                  of the PDB. The only differences from the standard 
                  are:
                  1. The residue name is 4 characters rather than 3 
                     (allowing LYSH, HISA, etc.).
                  2. The atom name starts one column later than the 
                     standard and is four columns wide encompasing the 
                     standard's `alternate' field. These two 
                     differences from the standard reflect the common
                     usage.
   V1.2a 28.06.90 Buffer size increased to 85 chars.
   V1.2b 15.02.91 Simply changed comment header to match new standard.
   V1.3  07.01.92 Corrected small bug in while() loop. Now ignores 
                  blank lines properly
   V1.4  11.05.92 Added check on EOF in while() loop and memset() of 
                  buffer. ANSIfied.
   V1.5  01.06.92 Documented for autodoc
   V1.7  01.10.92 Changed to use fgets()
   V1.6  19.06.92 Corrected use of stdlib
   V1.8  08.12.92 SAS/C V6 now defines atof() in stdlib
   V1.9  10.06.93 Returns TRUE or FALSE rather than exiting on failure
   V2.0  17.06.93 Rewritten to use fsscanf()
   V2.1  08.07.93 Modified to give ReadPDB() and ReadPDBAtoms()
   V2.2  09.07.93 Modified to return the PDB pointer rather than a BOOL.
                  There is now no need to initialise the structure first.
                  Rewrote allocation scheme.
   V2.3  17.03.94 Handles partial occupancies. If occupancies are not
                  1.0 or 0.0, the normal routine now reads only the 
                  highest occupancy atoms and corrects the atoms names 
                  to remove alternative labels. This behaviour can be 
                  overridden by calling one of the ...OccRank()
                  routines to read lower occupancy atoms. 
                  Sets natom to -1 if there was an error to distinguish 
                  from no atoms.
                  Handles atom names which start in column 13 rather
                  than column 14. This is allowed in the standard, but
                  very rare.
                  Added ReadPDBOccRank() & ReadPDBAtomsOccRank()
                  Sets gPDBPartialOcc flag.
   V2.4  06.04.94 With atom names which start in column 13, now checks
                  if the first character is a digit. If so, moves it
                  to the end of the atom name. Thus, 1HH1 becomes HH11
                  and 2HH1 becomes HH12.
   V2.5  04.10.94 Fixed partial occ when resnum changes as well as atom
                  name. Fixed bug when MAXPARTIAL exceeded.
   V2.6  03.11.94 Simply Corrected description. No code changes
   V2.7  06.03.95 Now reads just the first NMR model by default
                  doReadPDB() no longer static
                  Sets gPDBMultiNMR if ENDMDL records found.
   V2.8  13.01.97 Added check on return from fsscanf. Blank lines used
                  to result in duplication of the previous line since
                  fsscanf() does not reset the variables on receiving
                  a blank line. Also fixed in fsscanf().
   V2.9  25.02.98 Added transparent reading of gzipped PDB files if
                  GUNZIP_SUPPORT is defined
   V2.10 18.08.98 Added cast to popen() for SunOS
   V2.11 08.10.99 Initialised some variables
   V2.12 15.02.01 Added atnam_raw into PDB structure
   V2.13 30.05.02 Changed PDB field from 'junk' to 'record_type'
   V2.14 27.04.05 Fixed bug in atnam_raw for multiple occupancies
   V2.15 03.06.05 Added altpos field to PDB structure. The massaged atom
                  name no longer contains the alternate indicator and
                  atnam_raw has only the atom name with altpos having the
                  alternate indicator (as it should!)
   V2.16 14.10.05 Fixed a problem in StoreOccRankAtom() when a lower
                  occupancy atom has (erroneously) been set to occupancy
                  of zero and you want to pull out that atom
   V2.17 25.01.06 Added calls to RemoveAlternates()
   V2.18 03.02.06 Added prototypes for popen() and pclose()
   V2.19 05.06.07 Added support for Unix compress'd files
   V2.20 29.06.07 popen() and pclose() prototypes now skipped for MAC OSX
                  which defines them differently
   V2.21 17.03.09 popen() prototype skipped for Windows. By: CTP
   V2.22 21.12.11 doReadPDB() modified for cases where atoms are single
                  occupancy but occupancy is < 1.0

*************************************************************************/
/* Defines required for includes
*/
#define READPDB_MAIN

/************************************************************************/
/* Includes
*/
#include "port.h"    /* Required before stdio.h                         */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>

#include "SysDefs.h"
#include "MathType.h"
#include "pdb.h"
#include "macros.h"
#include "fsscanf.h"
#include "general.h"

#define MAXPARTIAL 8
#define SMALL      0.000001

/************************************************************************/
/* Prototypes
*/
static BOOL StoreOccRankAtom(int OccRank, PDB multi[MAXPARTIAL], 
                             int NPartial, PDB **ppdb, PDB **pp, 
                             int *natom);
#if !defined(__APPLE__) && !defined(MS_WINDOWS)
FILE *popen(char *, char *);
#endif
#ifndef __APPLE__
int  pclose(FILE *);
#endif

/************************************************************************/
/*>PDB *ReadPDB(FILE *fp, int *natom)
   ----------------------------------
   Input:   FILE     *fp      A pointer to type FILE in which the
                              .PDB file is stored.
   Output:  int      *natom   Number of atoms read. -1 if error.
   Returns: PDB      *pdb     A pointer to the first allocated item of
                              the PDB linked list

   Reads a PDB file into a PDB linked list

   08.07.93 Written as entry for doReadPDB()
   09.07.93 Modified to return pointer to PDB
   17.03.94 Modified to handle OccRank
   06.03.95 Added value for NMR model to read (1 = first)
   25.01.06 Added call to RemoveAlternates() - this deals with odd
            cases where alternate atom positions don't appear where
            they should!
   25.01.06 Added call to RemoveAlternates(). This deals with odd uses
            of multiple occupancies like 3pga and the instance where
            the alternates are all grouped at the end of the file.
*/
PDB *ReadPDB(FILE *fp,
             int  *natom)
{
   PDB *pdb;
   pdb = doReadPDB(fp, natom, TRUE, 1, 1);
   pdb = RemoveAlternates(pdb);
   return(pdb);
}

/************************************************************************/
/*>PDB *ReadPDBAll(FILE *fp, int *natom)
   -------------------------------------
   Input:   FILE     *fp      A pointer to type FILE in which the
                              .PDB file is stored.
   Output:  int      *natom   Number of atoms read. -1 if error.
   Returns: PDB      *pdb     A pointer to the first allocated item of
                              the PDB linked list

   Reads a PDB file into a PDB linked list. Reads all partial occupancy
   atoms. Reads both ATOM and HETATM records.

   04.10.94 Original    By: ACRM
   06.03.95 Added value for NMR model to read (0 = all)   
*/
PDB *ReadPDBAll(FILE *fp,
             int  *natom)
{
   return(doReadPDB(fp, natom, TRUE, 0, 0));
}

/************************************************************************/
/*>PDB *ReadPDBAtoms(FILE *fp, int *natom)
   ---------------------------------------
   Input:   FILE     *fp      A pointer to type FILE in which the
                              .PDB file is stored.
   Output:  int      *natom   Number of atoms read. -1 if error.
   Returns: PDB      *pdb     A pointer to the first allocated item of
                              the PDB linked list

   Reads a PDB file into a PDB linked list. Atoms only (no HETATM cards).

   08.07.93 Written as entry for doReadPDB()
   09.07.93 Modified to return pointer to PDB
   17.03.94 Modified to handle OccRank
   06.03.95 Added value for NMR model to read (1 = first)
   25.01.06 Added call to RemoveAlternates(). This deals with odd uses
            of multiple occupancies like 3pga and the instance where
            the alternates are all grouped at the end of the file.
*/
PDB *ReadPDBAtoms(FILE *fp,
                  int  *natom)
{
   PDB *pdb;
   pdb = doReadPDB(fp, natom, FALSE, 1, 1);
   pdb = RemoveAlternates(pdb);
   return(pdb);
}

/************************************************************************/
/*>PDB *ReadPDBOccRank(FILE *fp, int *natom, int OccRank)
   ------------------------------------------------------
   Input:   FILE     *fp      A pointer to type FILE in which the
                              .PDB file is stored.
            int      OccRank  Occupancy ranking (>=1)
   Output:  int      *natom   Number of atoms read. -1 if error.
   Returns: PDB      *pdb     A pointer to the first allocated item of
                              the PDB linked list

   Reads a PDB file into a PDB linked list selecting the OccRank'th
   highest occupancy atoms

   17.03.94 Original    By: ACRM
   06.03.95 Added value for NMR model to read (1 = first)
*/
PDB *ReadPDBOccRank(FILE *fp, int *natom, int OccRank)
{
   return(doReadPDB(fp, natom, TRUE, OccRank, 1));
}

/************************************************************************/
/*>PDB *ReadPDBAtomsOccRank(FILE *fp, int *natom, int OccRank)
   -----------------------------------------------------------
   Input:   FILE     *fp      A pointer to type FILE in which the
                              .PDB file is stored.
            int      OccRank  Occupancy ranking (>=1)
   Output:  int      *natom   Number of atoms read. -1 if error.
   Returns: PDB      *pdb     A pointer to the first allocated item of
                              the PDB linked list

   Reads a PDB file into a PDB linked list ignoring HETATM records
   and selecting the OccRank'th highest occupancy atoms

   17.03.94 Original    By: ACRM
   06.03.95 Added value for NMR model to read (1 = first)
*/
PDB *ReadPDBAtomsOccRank(FILE *fp, int *natom, int OccRank)
{
   return(doReadPDB(fp, natom, FALSE, OccRank, 1));
}

/************************************************************************/
/*>PDB *doReadPDB(FILE *fp, int *natom, BOOL AllAtoms, int OccRank,
                  int ModelNum)
   ----------------------------------------------------------------
   Input:   FILE     *fp      A pointer to type FILE in which the
                              .PDB file is stored.
            BOOL     AllAtoms TRUE:  ATOM & HETATM records
                              FALSE: ATOM records only
            int      OccRank  Occupancy ranking
            int      ModelNum NMR Model number (0 = all)
   Output:  int      *natom   Number of atoms read. -1 if error.
   Returns: PDB      *pdb     A pointer to the first allocated item of
                              the PDB linked list

   Reads a PDB file into a PDB linked list. The OccRank value indicates
   occupancy ranking to read for partial occupancy atoms.
   If any partial occupancy atoms are read the global flag 
   gPDBPartialOcc is set to TRUE.

   04.11.88 V1.0  Original
   07.02.89 V1.1  Ignore records which aren't ATOM or HETATM
   28.03.90 V1.2  Altered field widths to match PDB standard better
                  See notes above for deviations
   28.06.90 V1.2a Buffer size increased to 85 chars.
   15.02.91 V1.2b Changed comment header to match new standard.
   07.01.92 V1.3  Ignores blank lines properly
   11.05.92 V1.4  Check on EOF in while() loop, memset() buffer. 
                  ANSIed.
   01.06.92 V1.5  Documented for autodoc
   19.06.92 V1.6  Corrected use of stdlib
   01.10.92 V1.7  Changed to use fgets()
   10.06.93 V1.9  Returns 0 on failure rather than exiting
                  Replaced SIZE with sizeof(PDB) directly
   17.06.93 V2.0  Rewritten to use fsscanf()
   08.07.93 V2.1  Split from ReadPDB()
   09.07.93 V2.2  Modified to return pointer to PDB. Rewrote allocation
                  scheme.
   17.03.94 V2.3  Handles partial occupancies
                  Sets natom to -1 if there was an error to distinguish 
                  from no atoms.
                  Handles atom names which start in column 13 rather
                  than column 14. This is allowed in the standard, but
                  very rare.
                  Sets flag for partials.
   06.04.94 V2.4  Atom names starting in column 13 have their first
                  character moved to the end if it is a digit.
   03.10.94 V2.5  Check residue number as well as atom name when running
                  through alternative atoms for partial occupancy
                  Moved increment of NPartial, so only done if there
                  is space in the array. If OccRank is 0, all atoms are
                  read regardless of occupancy.
   06.03.95 V2.7  Added value for NMR model to read (0 = all)
                  No longer static. Sets gPDBMultiNMR if ENDMDL records
                  found.
   13.01.97 V2.8  Added check on return from fsscanf. Blank lines used
                  to result in duplication of the previous line since
                  fsscanf() does not reset the variables on receiving
                  a blank line. Also fixed in fsscanf().
   25.02.98 V2.9  Added code to read gzipped PDB files transparently
                  when GUNZIP_SUPPORT is defined
   17.08.98 V2.10 Added case to popen() for SunOS
   08.10.99 V2.11 Initialise CurIns and CurRes
   15.02.01 V2.12 Added atnam_raw
   27.04.05 V2.14 Added another atnam_raw for multiple occupancies
   03.06.05 V2.15 Added altpos
   14.10.05 V2.16 Modified detection of partial occupancy. handles
                  residues like 1zeh/B16 where a lower partial is
                  erroneously set to zero
   05.06.07 V2.19 Added support for Unix compress'd files
   21.12.11 V2.22 Modified for cases of single occupancy < 1.0
*/
PDB *doReadPDB(FILE *fpin,
               int  *natom,
               BOOL AllAtoms,
               int  OccRank,
               int  ModelNum)
{
   char     record_type[8],
            atnambuff[8],
            *atnam,
            atnam_raw[8],
            resnam[8],
            chain[4],
            insert[4],
            buffer[160],
            CurAtom[8],
            cmd[80],
            CurIns = ' ',
            altpos;
   int      atnum,
            resnum,
            CurRes = 0,
            NPartial,
            ModelCount = 1;
   FILE     *fp = fpin;
   double   x,y,z,
            occ,
            bval;
   PDB      *pdb  = NULL,
            *p,
            multi[MAXPARTIAL];   /* Temporary storage for partial occ   */

#ifdef GUNZIP_SUPPORT
   int      signature[3],
            i,
            ch;
#endif

   *natom         = 0;
   CurAtom[0]     = '\0';
   NPartial       = 0;
   gPDBPartialOcc = FALSE;
   gPDBMultiNMR   = FALSE;
   cmd[0]         = '\0';

#ifdef GUNZIP_SUPPORT
   /* See whether this is a gzipped file                                */
   for(i=0; i<3; i++)
      signature[i] = fgetc(fpin);
   for(i=2; i>=0; i--)
      ungetc(signature[i], fpin);
   if(((signature[0] == (int)0x1F) &&    /* gzip                        */
       (signature[1] == (int)0x8B) &&
       (signature[2] == (int)0x08)) ||
      ((signature[0] == (int)0x1F) &&    /* 05.06.07 compress           */
       (signature[1] == (int)0x9D) &&
       (signature[2] == (int)0x90)))
   {
      /* It is gzipped so we'll open gunzip as a pipe and send the data
         through that into a temporary file
      */
      sprintf(cmd,"gunzip >/tmp/readpdb_%d",(int)getpid());
      if((fp = (FILE *)popen(cmd,"w"))==NULL)
      {
         *natom = (-1);
         return(NULL);
      }
      while((ch=fgetc(fpin))!=EOF)
         fputc(ch, fp);
      pclose(fp);

      /* We now reopen the temporary file as our PDB input file         */
      sprintf(cmd,"/tmp/readpdb_%d",(int)getpid());
      if((fp = fopen(cmd,"r"))==NULL)
      {
         *natom = (-1);
         return(NULL);
      }
   }
#endif   

   while(fgets(buffer,159,fp))
   {
      if(ModelNum != 0)          /* We are interested in model numbers  */
      {
         if(!strncmp(buffer,"ENDMDL",6))
         {
            ModelCount++;
         }

         if(ModelCount < ModelNum)   /* Haven't reached the right model */
            continue;
         else if(ModelCount > ModelNum)    /* Gone past the right model */
            break;
      }

      if(!strncmp(buffer,"ENDMDL",6))
         gPDBMultiNMR   = TRUE;
      
      if(fsscanf(buffer,"%6s%5d%1x%5s%4s%1s%4d%1s%3x%8lf%8lf%8lf%6lf%6lf",
                 record_type,&atnum,atnambuff,resnam,chain,&resnum,insert,
                 &x,&y,&z,&occ,&bval) != EOF)
      {
         if((!strncmp(record_type,"ATOM  ",6)) || 
            (!strncmp(record_type,"HETATM",6) && AllAtoms))
         {
            /* Copy the raw atom name                                   */
            /* 03.06.05 Note: this reads the alternate atom position as 
               well as the atom name - changes in FixAtomName() now strip
               that
               We now copy only the first 4 characters into atnam_raw and
               put the 5th character into altpos
            */
            strncpy(atnam_raw, atnambuff, 4);
            atnam_raw[4] = '\0';
            altpos = atnambuff[4];

            /* Fix the atom name accounting for start in column 13 or 14*/
            atnam = FixAtomName(atnambuff, occ);
            
            /* Check for full occupancy. If occupancy is 0.0 assume that 
               it is actually fully occupied; the column just hasn't been
               filled in correctly
               
               04.10.94 Read all atoms if OccRank is 0

               14.10.05 Now takes an atom as full occupancy:
                           if occ==1.0
                           if occ==0.0 and altpos==' '
                           if OccRank==0
                        This fixes problems where a lower (partial)
                        occupancy has erroneously been set to zero
               21.12.11 Now only worries about partial occupancy if altpos
                        is a space. The first line of the if() statement
                        here would assume single occupancy if altpos was
                        a space and occupancy was zero:
                        if(((altpos == ' ') && (occ < (double)SMALL)) ||
                        - it now assumes single occupancy if altpos is a
                        space regardless of the actual occupancy. This
                        deals with cases like 1ap2 ZN A112 and 1ces ZN
                        A238 where these HETATMs are single occupancy
                        but with occupancy < 1.0
            */
            if((altpos == ' ') ||
               (occ > (double)0.999) || 
               (OccRank == 0))
            {
               /* Trim the atom name to 4 characters                    */
               atnam[4] = '\0';
               
               if(NPartial != 0)
               {
                  if(!StoreOccRankAtom(OccRank,multi,NPartial,&pdb,&p,
                                       natom))
                  {
                     if(pdb != NULL) FREELIST(pdb, PDB);
                     *natom = (-1);
                     if(cmd[0]) unlink(cmd);
                     return(NULL);
                  }
                  
                  /* Set partial occupancy counter to 0                 */
                  NPartial = 0;
               }
               
               /* Allocate space in the linked list                     */
               if(pdb == NULL)
               {
                  INIT(pdb, PDB);
                  p = pdb;
               }
               else
               {
                  ALLOCNEXT(p, PDB);
               }
               
               /* Failed to allocate space; free up list so far & return*/
               if(p==NULL)
               {
                  if(pdb != NULL) FREELIST(pdb, PDB);
                  *natom = (-1);
                  if(cmd[0]) unlink(cmd);
                  return(NULL);
               }
               
               /* Increment the number of atoms                         */
               (*natom)++;
               
               /* Store the information read                            */
               p->atnum  = atnum;
               p->resnum = resnum;
               p->x      = (REAL)x;
               p->y      = (REAL)y;
               p->z      = (REAL)z;
               p->occ    = (REAL)occ;
               p->bval   = (REAL)bval;
               p->altpos = altpos;    /* 03.06.05 Added this one        */
               p->next   = NULL;
               strcpy(p->record_type, record_type);
               strcpy(p->atnam,       atnam);
               strcpy(p->atnam_raw,   atnam_raw);
               strcpy(p->resnam,      resnam);
               strcpy(p->chain,       chain);
               strcpy(p->insert,      insert);
            }
            else   /* Partial occupancy                                 */
            {
               /* Set flag to say we've got a partial occupancy atom    */
               gPDBPartialOcc = TRUE;
               
               /* First in a group, store atom name                     */
               if(NPartial == 0)
               {
                  CurIns = insert[0];
                  CurRes = resnum;
                  strncpy(CurAtom,atnam,8);
               }
               
               if(strncmp(CurAtom,atnam,strlen(CurAtom)-1) || 
                  resnum != CurRes || 
                  CurIns != insert[0])
               {
                  /* Atom name has changed 
                     Select and store the OccRank highest occupancy atom
                     */
                  if(!StoreOccRankAtom(OccRank,multi,NPartial,&pdb,&p,
                                       natom))
                  {
                     if(pdb != NULL) FREELIST(pdb, PDB);
                     *natom = (-1);
                     if(cmd[0]) unlink(cmd);
                     return(NULL);
                  }
                  
                  /* Reset the partial atom counter                     */
                  NPartial = 0;
                  strncpy(CurAtom,atnam,8);
                  CurRes = resnum;
                  CurIns = insert[0];
               }
               
               if(NPartial < MAXPARTIAL)
               {
                  /* Store the partial atom data                        */
                  multi[NPartial].atnum  = atnum;
                  multi[NPartial].resnum = resnum;
                  multi[NPartial].x      = (REAL)x;
                  multi[NPartial].y      = (REAL)y;
                  multi[NPartial].z      = (REAL)z;
                  multi[NPartial].occ    = (REAL)occ;
                  multi[NPartial].bval   = (REAL)bval;
                  multi[NPartial].next   = NULL;
                  strcpy(multi[NPartial].record_type, record_type);
                  strcpy(multi[NPartial].atnam,       atnam);
                  /* 27.04.05 - added this line                         */
                  strcpy(multi[NPartial].atnam_raw,   atnam_raw);
                  strcpy(multi[NPartial].resnam,      resnam);
                  strcpy(multi[NPartial].chain,       chain);
                  strcpy(multi[NPartial].insert,      insert);
                  /* 03.06.05 - added this line                         */
                  multi[NPartial].altpos = altpos;
                  
                  NPartial++;
               }
            }
         }
      }
   }

   if(NPartial != 0)
   {
      if(!StoreOccRankAtom(OccRank,multi,NPartial,&pdb,&p,natom))
      {
         if(pdb != NULL) FREELIST(pdb, PDB);
         *natom = (-1);
         if(cmd[0]) unlink(cmd);
         return(NULL);
      }
   }

   if(cmd[0]) unlink(cmd);

   /* Return pointer to start of linked list                            */
   return(pdb);
}

/************************************************************************/
/*>static BOOL StoreOccRankAtom(int OccRank, PDB multi[MAXPARTIAL], 
                                int NPartial, PDB **ppdb, PDB **pp, 
                                int *natom)
   ----------------------------------------------------------------
   Input:   int  OccRank     Occupancy ranking required (>=1)
            PDB  multi[]     Array of PDB records for alternative atom
                             positions
            int  NPartial    Number of items in multi array
   I/O:     PDB  **ppdb      Start of PDB linked list (or NULL)
            PDB  **pp        Current position in PDB linked list (or NULL)
            int  *natom      Number of atoms read
   Returns: BOOL             Memory allocation success

   Takes an array of PDB records which represent alternative atom 
   positions for an atom. Select the OccRank'th highest occupancy and
   add this one into the PDB linked list.

   To be called by doReadPDB().

   17.03.94 Original    By: ACRM
   08.10.99 Initialise IMaxOcc and MaxOcc
   27.04.05 Added atnam_raw
   03.06.05 Added altpos
   14.10.05 Modified the flag value from 0.0 to -1.0 so that erroneous
            lower occupancies of 0.0 are read properly and written back
            with their occupancy (0.0) rather than the next higher
            occupancy. Handles residues like 1zeh/B16
*/
static BOOL StoreOccRankAtom(int OccRank, PDB multi[MAXPARTIAL], 
                             int NPartial, PDB **ppdb, PDB **pp, 
                             int *natom)
{
   int  i,
        j,
        IMaxOcc = 0;
   REAL MaxOcc  = (REAL)0.0,
        LastOcc = (REAL)0.0;
   
   if(OccRank < 1) OccRank = 1;
   
   for(i=0; i<OccRank; i++)
   {
      MaxOcc  = (REAL)0.0;
      IMaxOcc = 0;
      
      for(j=0; j<NPartial; j++)
      {
         if(multi[j].occ >= MaxOcc)
         {
            MaxOcc  = multi[j].occ;
            IMaxOcc = j;
         }
      }
      /* 14.10.05 Changed flag value to -1 so that erroneous occupancies
         of zero are treated properly
      */
      multi[IMaxOcc].occ = (REAL)-1.0;

      /* 14.10.05 Changed flag value to -1 so that erroneous occupancies
         of zero are treated properly
      */
      if(MaxOcc < (REAL)0.0) break;
      LastOcc = MaxOcc;
   }

   /* If we ran out of rankings, take the last one to be found          */
   /* 14.10.05 Changed flag value to -1 so that erroneous occupancies
      of zero are treated properly
   */
   if(MaxOcc < (REAL)0.0)
      MaxOcc = LastOcc;

   /* Store this atom
      Allocate space in the linked list
   */
   if(*ppdb == NULL)
   {
      INIT((*ppdb), PDB);
      *pp = *ppdb;
   }
   else
   {
      ALLOCNEXT(*pp, PDB);
   }
            
   /* Failed to allocate space; error return.                           */
   if(*pp==NULL)
      return(FALSE);
               
   /* Increment the number of atoms                                     */
   (*natom)++;
               
   /* Store the information read                                        */
   (*pp)->atnum  = multi[IMaxOcc].atnum;
   (*pp)->resnum = multi[IMaxOcc].resnum;
   (*pp)->x      = multi[IMaxOcc].x;
   (*pp)->y      = multi[IMaxOcc].y;
   (*pp)->z      = multi[IMaxOcc].z;
   (*pp)->occ    = MaxOcc;
   (*pp)->bval   = multi[IMaxOcc].bval;
   (*pp)->next   = NULL;
   /* 03.06.05 Added this line                                          */
   (*pp)->altpos = multi[IMaxOcc].altpos;
   strcpy((*pp)->record_type, multi[IMaxOcc].record_type);
   strcpy((*pp)->atnam,       multi[IMaxOcc].atnam);
   /* 27.04.05 Added this line                                          */
   strcpy((*pp)->atnam_raw,   multi[IMaxOcc].atnam_raw);
   strcpy((*pp)->resnam,      multi[IMaxOcc].resnam);
   strcpy((*pp)->chain,       multi[IMaxOcc].chain);
   strcpy((*pp)->insert,      multi[IMaxOcc].insert);

   /* Patch the atom name to remove the alternate letter                */
   if(strlen((*pp)->atnam) > 4)
      ((*pp)->atnam)[4] = '\0';
   else
      ((*pp)->atnam)[3] = ' ';

   return(TRUE);
}

/************************************************************************/
/*>char *FixAtomName(char *name, REAL occup)
   -----------------------------------------
   Input:   char  *name     Atom name read from file
            REAL  occup     Occupancy to allow fixing of partial occupancy
                            atom names
   Returns: char  *         Fixed atom name (pointer into name)

   Fixes an atom name by removing leading spaces, or moving a leading
   digit to the end of the string. Used by doReadPDB()

   06.04.94 Original    By: ACRM
   01.03.01 No longer static
   03.06.05 The name passed in has always contained the column which is
            officially the alternate atom position indicator, but is 
            used by some programs as part of the atom name. Thus the
            properly constructed variable coming into the routine should
            be something like '1HG1 ' or '1HG1A' for an alternate atom
            position. However some programs use ' HG11'. Therefore we
            now check for a character in the last position and replace
            it with a space if there is a space in the preceeding
            position (e.g. ' CA A' -> ' CA  ') or if there is a 
            character in the first position (e.g. '1HG1A' -> '1HG1 ')
            or if the occupancy is not zero/one

            NOTE!!! To support this, the routine now has a second 
            parameter: REAL occup
*/
char *FixAtomName(char *name, REAL occup)
{
   char *newname;
   int  len;

   /* Default behaviour, just return the input string                   */
   newname = name;

   if(name[0] == ' ')       /* Name starts in column 14                 */
   {
      /* remove leading spaces                                          */
      KILLLEADSPACES(newname,name);
      /* 03.06.05 If the last-but-one position is a space, force the last
         position (the alternate atom indicator) to be a space
      */
      if(newname[2] == ' ')
      {
         newname[3] = ' ';
      }
   }
   else                     /* Name starts in column 13                 */
   {
      /* 03.06.05 The last character is the alternate atom indicator, 
         so force it to be a space
      */
      name[4] = ' ';
      
      /* If the first character is a digit, move it to the end          */
      if(isdigit(name[0]))
      {
         if((len = chindex(name,' ')) == (-1))
         {
            /* We didn't find a space in the name, so add the character
               onto the end of the string and re-terminate
            */
            len         = strlen(name);
            newname     = name+1;
            name[len]   = name[0];
            name[len+1] = '\0';
         }
         else
         {
            /* We did find a space in the name, so put the first
               character there
            */
            newname     = name+1;
            name[len]   = name[0];
         }
      }
   }
   return(newname);
}

/************************************************************************/
/*>PDB *RemoveAlternates(PDB *pdb)
   -------------------------------
   I/O:       PDB     *pdb       PDB 
   Returns:   PDB *              Ammended linked list (in case start has
                                 changed)

   Remove alternate atoms - we keep only the highest occupancy or the 
   first if there are more than one the same.

   25.01.05 Original based on code written for Inpharmatica   By: ACRM
*/
PDB *RemoveAlternates(PDB *pdb)
{
   PDB   *p, 
         *q, 
         *r, 
         *s, 
         *s_prev,
         *r_prev,
         *a_prev,
         *next,
         *alts[MAXPARTIAL];
   int   i, 
         altCount, 
         highest;
   
   
   /* Step through residues                                             */
   r_prev=NULL;
   for(p=pdb; p!=NULL; p=q)
   {
      q=FindNextResidue(p);
      
      /* Step through atoms                                             */
      for(r=p; r!=q; NEXT(r))
      {
         if(r->altpos != ' ')
         {
#ifdef DEBUG
            fprintf(stderr,"\n\nAlt pos found for record:\n");
            WritePDBRecord(stderr, r);
#endif
            /* We have an alternate, store it and search for the other 
               ones  
            */
            altCount=0;
            alts[altCount++] = r;
            /* Search through this residue for the alternates.
               This will work for 99.9% of files where the alternates are
               with the main atoms
            */
            for(s=r->next; s!=q; NEXT(s))
            {
               if(!strcmp(s->atnam_raw, alts[0]->atnam_raw))
               {
                  if(altCount < MAXPARTIAL)
                  {
                     alts[altCount++] = s;
#ifdef DEBUG
                     fprintf(stderr,"Partner atom found in res:\n");
                     WritePDBRecord(stderr, s);
#endif
                  }
                  else
                  {
                     fprintf(stderr,"Warning==> More than %d alternative \
conformations in\n", MAXPARTIAL);
                     fprintf(stderr,"           residue %c%d%c atom %s. \
Increase MAXPARTIAL in ReadPDB.c\n", s->chain[0],
                                     s->resnum,
                                     s->insert[0],
                                     s->atnam);
                  }
               }
            }
            /* If we didn't find the alternates within the residue, then
               we search the rest of the records.
               This covers the known entry where the alternates are shoved
               on the end instead! 
            */
            if(altCount<2)
            {
#ifdef DEBUG
               fprintf(stderr,"No partner found in residue\n");
#endif

               s_prev = NULL;
               for(s=q; s!=NULL; NEXT(s))
               {
                  if((s->resnum    == alts[0]->resnum) &&
                     (s->insert[0] == alts[0]->insert[0]) &&
                     (s->chain[0]  == alts[0]->chain[0]) &&
                     !strcmp(s->atnam_raw, alts[0]->atnam_raw))
                  {
                     if(altCount < MAXPARTIAL)
                     {
                        alts[altCount++] = s;
#ifdef DEBUG
                        fprintf(stderr,"Partner found outside \
residue:\n");
                        WritePDBRecord(stderr, s);
#endif
                     }
                     else
                     {
                        fprintf(stderr,"Warning==> More than %d \
alternative conformations in\n", MAXPARTIAL);
                        fprintf(stderr,"           residue %c%d%c atom \
%s. Increase MAXPARTIAL in ReadPDB.c\n", s->chain[0],
                                         s->resnum,
                                         s->insert[0],
                                         s->atnam);

                        /* Move this record to the correct position in the
                           linked list

                           First unlink s from its old position
                        */
                        if(s_prev != NULL)
                           s_prev->next = s->next;
                     
                        /* Now link it back in where it should be       */
                        next = r->next;
                        r->next = s;
                        s->next = next;
                     }
                  }
                  s_prev = s;
               }
            }

            if(altCount < 2)
            {
#ifdef DEBUG
               fprintf(stderr,"No alternates found. Resetting ALT \
flag\n\n");
#endif
               alts[0]->altpos = ' ';
               
            }
            else
            {
               /* Find the highest occupancy, defaulting to the first   */
               highest = 0;
               for(i=0; i<altCount; i++)
               {
                  if(alts[i]->occ > alts[highest]->occ)
                     highest = i;
               }
               
               /* Delete the unwanted alternates                        */
               for(i=0; i<altCount; i++)
               {
                  if(i==highest) /* For the highest remove the ALT flag */
                  {
#ifdef DEBUG
                     fprintf(stderr,"Highest occupancy selected:\n");
                     WritePDBRecord(stderr, alts[i]);
#endif
                     alts[i]->altpos = ' ';
                  }
                  else
                  {
                     /* If we are deleting the current record pointer, 
                        then we need to update it
                     */
                     if(alts[i] == r)
                     {
#ifdef DEBUG
                        fprintf(stderr,"Deleting current record \
pointer\n");
#endif
                        
                        if(r_prev == NULL)
                        {
                           r_prev = r;
                           NEXT(r);
                           /* We are deleting the head of the list so we 
                              must update the main list pointer
                           */
                           pdb = r;
                        }
                        else
                        {
                           r = r_prev;
                           FINDPREV(r_prev, pdb, r);
                        }
                     }
                     
                     /* Delete the alternate we don't need              */
#ifdef DEBUG
                     fprintf(stderr,"Deleting Alt pos record:\n");
                     WritePDBRecord(stderr, alts[i]);
#endif
                     
                     FINDPREV(a_prev, pdb, alts[i]);
                     if(a_prev != NULL)
                        a_prev->next = alts[i]->next;
                     free(alts[i]);
                     
                  }  /* Not the highest, so we delete it                */
               }  /* Stepping through the alternates                    */
            }
            
         }  /* We have an alternate                                     */
         r_prev = r;
      }  /* Stepping through the atoms of this residue                  */
   }  /* Stepping through the residues                                  */
   return(pdb);
}

