/*************************************************************************

   Program:    
   File:       SplitSeq.c
   
   Version:    V1.10
   Date:       02.10.00
   Function:   
   
   Copyright:  (c) SciTech Software 1993-2000
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
   V1.0  29.09.92 Original
   V1.1  07.06.93 Corrected allocation
   V1.2  18.06.93 Handles multi-chains and skips NTER and CTER residues.
                  Added SplitSeq()
   V1.3  09.07.93 SplitSeq() cleans up properly if allocation failed
   V1.4  11.05.94 Added TrueSeqLen()
   V1.5  13.05.94 Fixed bug in PDB2Seq().
                  Added KnownSeqLen().
   V1.6  07.09.94 Fixed allocation bug in SplitSeq()
   V1.7  19.07.95 Added check for ATOM records
   V1.8  24.01.96 Fixed bug when no ATOM records in linked list
                  Returns a blank string
   V1.9  26.08.97 Renamed DoPDB2Seq() with handling of Asx/Glx and
                  protein-only. Added macros to recreate the
                  old PDB2Seq() interface and similar new calls
   V1.10 02.10.00 Added NoX option

*************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include <string.h>

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
/*>int SplitSeq(char *LinearSeq, char **seqs)
   ------------------------------------------
   Input:   char  *LinearSeq   Array containing sequence with chains
                               terminated by *'s
   Output:  char  **seqs       Allocated set of character arrays 
                               containing one chain per array
   Returns: int                Number of chains found

   Splits a sequence stored as a linear array with each chain separated
   by a * into an array of sequences. Returns the number of chains
   found.
   
   18.06.93 Original    By: ACRM
   09.07.93 Cleans up properly of allocation failed
   07.09.94 Sequence space was being allocated one too small
*/
int SplitSeq(char *LinearSeq, char **seqs)
{
   char  *ptr,
         *star;
   int   seqlen   = strlen(LinearSeq),
         NSeq     = 0;
   
   ptr = LinearSeq;
   
   while(ptr-LinearSeq < seqlen)
   {
      star = strchr(ptr,'*');
      if(star != NULL)
      {
         *star = '\0';
         seqs[NSeq] = (char *)malloc((1+strlen(ptr)) * sizeof(char));
         if(seqs[NSeq] == NULL)
         {
            int i;
            for(i=0; i<=NSeq; i++) if(seqs[i] != NULL) free(seqs[i]);
            return(0);
         }
         strcpy(seqs[NSeq],ptr);
         NSeq++;
         ptr = star+1;
      }
      else
      {
         if(ptr < LinearSeq+seqlen)
         {
            seqs[NSeq] = (char *)malloc((1+strlen(ptr)) * sizeof(char));
            if(seqs[NSeq] == NULL)
            {
               int i;
               for(i=0; i<=NSeq; i++) if(seqs[i] != NULL) free(seqs[i]);
               return(0);
            }
            strcpy(seqs[NSeq],ptr);
            NSeq++;
         }
         break;
      }
   }
   
   return(NSeq);
}

