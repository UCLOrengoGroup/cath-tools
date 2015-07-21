/*************************************************************************

   Program:    
   File:       ReadWholePDB.c
   
   Version:    V1.3
   Date:       17.03.09
   Function:   
   
   Copyright:  (c) Dr. Andrew C. R. Martin, University of Reading, 2002
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
   V1.0  30.05.02 Original
   V1.1  12.06.08 CTP Added include for port.h
   V1.2  13.06.08 popen() and pclose() prototypes skipped for Mac OS X.
   V1.3  17.03.09 popen() prototype skipped for Windows. By: CTP

*************************************************************************/
/* Includes
*/
#include "port.h"    /* Required before stdio.h                         */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "macros.h"
#include "general.h"
#include "pdb.h"


/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160


/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static WHOLEPDB *doReadWholePDB(FILE *fpin, BOOL atomsonly);

#if !defined(__APPLE__) && !defined(MS_WINDOWS)
FILE *popen(char *, char *);
#endif
#ifndef __APPLE__
int  pclose(FILE *);
#endif

/************************************************************************/
/*>void FreeWholePDB(WHOLEPDB *wpdb)
   ---------------------------------
   Input:     WHOLEPDB    *wpdb    WHOLEPDB structure to be freed

   Frees the header, trailer and atom content from a WHOLEPDB structure

   30.05.02  Original   By: ACRM
*/
void FreeWholePDB(WHOLEPDB *wpdb)
{
   FreeStringList(wpdb->header);
   FreeStringList(wpdb->trailer);
   FREELIST(wpdb->pdb, PDB);
   free(wpdb);
}

/************************************************************************/
/*>void WriteWholePDB(FILE *fp, WHOLEPDB *wpdb)
   --------------------------------------------
   Input:     FILE       *fp        File pointer
              WHOLEPDB   *wpdb      Whole PDB structure pointer

   Writes a PDB file including header and trailer information

   30.05.02  Original   By: ACRM
*/
void WriteWholePDB(FILE *fp, WHOLEPDB *wpdb)
{
   WriteWholePDBHeader(fp, wpdb);
   WritePDB(fp, wpdb->pdb);
   WriteWholePDBTrailer(fp, wpdb);
}


/************************************************************************/
/*>void WriteWholePDBHeader(FILE *fp, WHOLEPDB *wpdb)
   --------------------------------------------------
   Input:     FILE       *fp        File pointer
              WHOLEPDB   *wpdb      Whole PDB structure pointer

   Writes the header of a PDB file 

   30.05.02  Original   By: ACRM
*/
void WriteWholePDBHeader(FILE *fp, WHOLEPDB *wpdb)
{
   STRINGLIST *s;
   
   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      fputs(s->string, fp);
   }
}


/************************************************************************/
/*>void WriteWholePDBTrailer(FILE *fp, WHOLEPDB *wpdb)
   ---------------------------------------------------
   Input:     FILE       *fp        File pointer
              WHOLEPDB   *wpdb      Whole PDB structure pointer

   Writes the trailer of a PDB file 

   30.05.02  Original   By: ACRM
*/
void WriteWholePDBTrailer(FILE *fp, WHOLEPDB *wpdb)
{
   STRINGLIST *s;
   
   for(s=wpdb->trailer; s!=NULL; NEXT(s))
   {
      fputs(s->string, fp);
   }
}


/************************************************************************/
/*>WHOLEPDB *ReadWholePDB(FILE *fpin)
   ----------------------------------
   Input:     FILE      *fpin     File pointer
   Returns:   WHOLEPDB  *         Whole PDB structure containing linked
                                  list to PDB coordinate data

   Reads a PDB file, storing the header and trailer information as
   well as the coordinate data. Can read gzipped files as well as
   uncompressed files.

   Coordinate data is accessed as linked list of type PDB as follows:
   
   WHOLEPDB *wpdb;
   PDB      *p;
   wpdb = ReadWholePDB(fp);
   for(p=wpdb->pdb; p!=NULL; p=p->next)
   {
      ... Do something with p ...
   }

   07.03.07 Made into a wrapper to doReadWholePDB()
*/
WHOLEPDB *ReadWholePDB(FILE *fpin)
{
   return(doReadWholePDB(fpin, FALSE));
}

/************************************************************************/
/*>WHOLEPDB *ReadWholePDBAtoms(FILE *fpin)
   ---------------------------------------
   Input:     FILE      *fpin     File pointer
   Returns:   WHOLEPDB  *         Whole PDB structure containing linked
                                  list to PDB coordinate data

   Reads a PDB file, storing the header and trailer information as
   well as the coordinate data. Can read gzipped files as well as
   uncompressed files.

   Coordinate data is accessed as linked list of type PDB as follows:
   
   WHOLEPDB *wpdb;
   PDB      *p;
   wpdb = ReadWholePDB(fp);
   for(p=wpdb->pdb; p!=NULL; p=p->next)
   {
      ... Do something with p ...
   }

   07.03.07 Made into a wrapper to doReadWholePDB()
*/
WHOLEPDB *ReadWholePDBAtoms(FILE *fpin)
{
   return(doReadWholePDB(fpin, TRUE));
}


/************************************************************************/
/*>static WHOLEPDB *ReadWholePDB(FILE *fpin)
   -----------------------------------------
   Input:     FILE      *fpin     File pointer
   Returns:   WHOLEPDB  *         Whole PDB structure containing linked
                                  list to PDB coordinate data

   Reads a PDB file, storing the header and trailer information as
   well as the coordinate data. Can read gzipped files as well as
   uncompressed files.

   Coordinate data is accessed as linked list of type PDB as follows:
   
   WHOLEPDB *wpdb;
   PDB      *p;
   wpdb = ReadWholePDB(fp);
   for(p=wpdb->pdb; p!=NULL; p=p->next)
   {
      ... Do something with p ...
   }

   30.05.02 Original   By: ACRM
   07.03.07 Made into a doXXX routine to add a atomsonly parameter
   05.06.07 Added support for Unix compress'd files

   TODO FIXME!!!!! Move all this into doReadPDB so that we don't worry 
   about rewinding any more
*/
static WHOLEPDB *doReadWholePDB(FILE *fpin, BOOL atomsonly)
{
   WHOLEPDB *wpdb;
   char     buffer[MAXBUFF];
   FILE     *fp = fpin;
   
#ifdef GUNZIP_SUPPORT
   int      signature[3],
            i,
            ch;
   char     cmd[80];
#endif

   if((wpdb=(WHOLEPDB *)malloc(sizeof(WHOLEPDB)))==NULL)
      return(NULL);

   wpdb->pdb     = NULL;
   wpdb->header  = NULL;
   wpdb->trailer = NULL;
   
#ifdef GUNZIP_SUPPORT
   cmd[0] = '\0';
   
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
      cmd[0] = '\0';
      sprintf(cmd,"gunzip >/tmp/readpdb_%d",(int)getpid());
      if((fp = (FILE *)popen(cmd,"w"))==NULL)
      {
         wpdb->natoms = (-1);
         return(NULL);
      }
      while((ch=fgetc(fpin))!=EOF)
         fputc(ch, fp);
      pclose(fp);

      /* We now reopen the temporary file as our PDB input file         */
      sprintf(cmd,"/tmp/readpdb_%d",(int)getpid());
      if((fp = fopen(cmd,"r"))==NULL)
      {
         wpdb->natoms = (-1);
         return(NULL);
      }
   }
#endif   

   /* Read the header from the PDB file                                 */
   while(fgets(buffer,MAXBUFF,fp))
   {
      if(!strncmp(buffer, "ATOM  ", 6) ||
         !strncmp(buffer, "HETATM", 6) ||
         !strncmp(buffer, "MODEL ", 6))
      {
         break;
      }
      if((wpdb->header = StoreString(wpdb->header, buffer))==NULL)
         return(NULL);
   }
   
   /* Read the coordinates                                              */
   rewind(fp);
   if(atomsonly)
   {
      wpdb->pdb = ReadPDBAtoms(fp, &(wpdb->natoms));
   }
   else
   {
      wpdb->pdb = ReadPDB(fp, &(wpdb->natoms));
   }

   /* Read the trailer                                                  */
   rewind(fp);
   while(fgets(buffer,MAXBUFF,fp))
   {
      if(!strncmp(buffer, "CONECT", 6) ||
         !strncmp(buffer, "MASTER", 6) ||
         !strncmp(buffer, "END   ", 6))
      {
         wpdb->trailer = StoreString(wpdb->trailer, buffer);
      }
   }
   
   return(wpdb);
}

