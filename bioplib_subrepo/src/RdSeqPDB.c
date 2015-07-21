/*************************************************************************

   Program:    
   File:       RdSeqPDB.c
   
   Version:    V1.0
   Date:       14.10.96
   Function:   Read sequence from SEQRES records in a PDB file
   
   Copyright:  (c) SciTech Software 1996
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      andrew@stagleys.demon.co.uk
               
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

*************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include "general.h"
#include "seq.h"
#include "macros.h"
#include "fsscanf.h"

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
char **ReadSeqresPDB(FILE *fp, int *nchains);
static STRINGLIST *RdSeqRes(FILE *fp);


/************************************************************************/
/*>char **ReadSeqresPDB(FILE *fp, int *nchains)
   --------------------------------------------
   Input:   FILE   *fp       PDB file pointer
   Output:  int    *nchain   Number of chains found
   Returns: char   **        Array of sequence strings

   Reads the sequence from the SEQRES records of a PDB file. Creates
   an array of malloc()'d character arrays in which the sequence is
   stored. Can therefore cope with any size of sequence information
   from the PDB file.

   This is not normally recommended to get the sequence for a PDB file
   this way, but is useful to detect discrepancies compared with the
   sequence described by the ATOM records.

   14.10.96 Original   By: ACRM
*/
char **ReadSeqresPDB(FILE *fp, int *nchains)
{
   STRINGLIST *seqres = NULL, 
              *s;
   char       currchain,
              chain,
              **seqs,
              res[13][8];
   int        chainnum = 0,
              nres     = 0,
              i;

   *nchains = 0;
   
   /* First read the SEQRES records into a linked list                  */
   if((seqres = RdSeqRes(fp))==NULL)
      return(NULL);

   /* FIRST PASS: See how many chains there are                         */
   currchain = seqres->string[11];
   *nchains  = 1;
   for(s=seqres; s!=NULL; NEXT(s))
   {
      chain = s->string[11];
      if(chain != currchain)
      {
         currchain = chain;
         (*nchains)++;
      }
   }

   /* Allocate an array of character pointers to store this number of
      strings
   */
   if((seqs=(char **)malloc((*nchains) * sizeof(char *)))==NULL)
   {
      FREELIST(seqres, STRINGLIST);
      return(NULL);
   }

   /* SECOND PASS: Allocate space to store each chain                   */
   chainnum  = 0;
   currchain = '\0';
   for(s=seqres; s!=NULL; NEXT(s))
   {
      fsscanf(s->string,"%11x%c%5d",&chain,&nres);
      if(chain != currchain)
      {
         currchain = chain;
         if((seqs[chainnum]=(char *)malloc((nres+1)*sizeof(char))) 
            == NULL)
         {
            FREELIST(seqres, STRINGLIST);
            return(NULL);
         }
         chainnum++;
      }
   }

   /* THIRD PASS: Store the sequence                                    */
   chainnum  = 0;
   nres      = 0;
   currchain = seqres->string[11];
   for(s=seqres; s!=NULL; NEXT(s))
   {
      fsscanf(s->string,"%11x%c%7x%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s",
              &chain,res[0],res[1],res[2],res[3],res[4],res[5],res[6],
              res[7],res[8],res[9],res[10],res[11],res[12]);
      if(chain != currchain)
      {
         /* Start of new chain, terminate last one                      */
         seqs[chainnum][nres] = '\0';
         currchain = chain;
         nres      = 0;
         chainnum++;
      }
      
      /* Store these sequence data                                      */
      for(i=0; i<13; i++)
      {
         /* Break out if not all positions were filled in               */
         if(res[i][0] == ' ')
            break;
         seqs[chainnum][nres++] = throne(res[i]);
      }
   }
   /* Terminate last chain                                              */
   seqs[chainnum][nres] = '\0';

   FREELIST(seqres, STRINGLIST);
   
   return(seqs);
}


/************************************************************************/
/*>static STRINGLIST *RdSeqRes(FILE *fp)
   -------------------------------------
   Input:   FILE        *fp      PDB File pointer
   Returns: STRINGLIST  *        Linked list of SEQRES records

   Used by ReadSeqresPDB() to read the SEQRES records into a linked list.

   14.10.96 Original   By: ACRM
*/
static STRINGLIST *RdSeqRes(FILE *fp)
{
   static STRINGLIST *seqres = NULL;
   char              buffer[MAXBUFF];
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      if(!strncmp(buffer,"SEQRES",6))
      {
         if((seqres = StoreString(seqres, buffer)) == NULL)
         {
            FREELIST(seqres, STRINGLIST);
            return(NULL);
         }
      }
   }
   
   return(seqres);
}

