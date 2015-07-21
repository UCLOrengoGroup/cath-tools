/*************************************************************************

   Program:    
   File:       ReadPIR.c
   
   Version:    V2.7R
   Date:       06.02.96
   Function:   Read a PIR sequence file
   
   Copyright:  (c) SciTech Software 1991-6
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

   int ReadPIR(FILE *fp, BOOL DoInsert, char **seqs, int maxchain, 
               SEQINFO *seqinfo, BOOL *punct, BOOL *error)
   ---------------------------------------------------------------
   This version attempts to read any PIR file following the PIR 
   specifications. It also accepts a few non-standard features:
   lower case sequence, no star at end of last chain, dashes in the
   sequence to indicate insertions.

   See also:
   int SimpleReadPIR(FILE *fp, int maxres, char **seqs)
   int ReadRawPIR(FILE *fp, BOOL DoInsert, char **seqs, int maxchain, 
                  SEQINFO *seqinfo, BOOL *punct, BOOL *error)

**************************************************************************

   Revision History:
   =================
   V1.0  01.06.92 Original
   V2.0  08.03.94 Changed name of ReadPIR() to ReadSimplePIR()
                  Added new ReadPIR().
   V2.1  18.03.94 getc() -> fgetc()
   V2.2  11.05.94 Changes to ReadPIR() for better compatibility with
                  PIR V38.0 and V39.0
   V2.3  28.02.95 Added ReadRawPIR()
   V2.4  13.03.95 Fixed bug in reading text lines in ReadRawPIR()
   V2.5  26.07.95 Removed unused variables
   V2.6  30.10.95 Cosmetic
   V2.7  06.02.96 Removes trailing spaces from comment line

*************************************************************************/
/* Includes
*/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "SysDefs.h"
#include "macros.h"
#include "seq.h"

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
/*>int ReadPIR(FILE *fp, BOOL DoInsert, char **seqs, int maxchain, 
               SEQINFO *seqinfo, BOOL *punct, BOOL *error)
   ------------------------------------------------------------------
   Input:   FILE    *fp      File pointer
            BOOL    DoInsert TRUE Read - characters into the sequence
                             FALSE Skip - characters
            int     maxchain Max number of chains to read. This is the
                             dimension of the seqs array.
                             N.B. THIS SHOULD BE AT LEAST 1 MORE THAN
                             THE EXPECTED MAXIMUM NUMBER OF SEQUENCES
   Output:  char    **seqs   Array of character pointers which will
                             be filled in with sequence information.
                             Memory will be allocated for any sequence
                             length.
            SEQINFO *seqinfo This structure will be filled in with
                             extra information about the sequence.
                             Header & title information and details
                             of any punctuation.
            BOOL    *punct   TRUE if any punctuation found.
            BOOL    *error   TRUE if an error occured (e.g. memory
                             allocation)
   Returns: int              Number of chains in this sequence.
                             0 if file ended, or no valid sequence
                             entries found.

   This is an all-singing, all-dancing PIR reader which should handle
   all legal PIR files and some (slightly) incorrect ones. The only
   requirements of the code are that the PIR file should have 2 title
   lines per entry, the first line starting with a > sign.

   The routine will handle multiple sequence files. Successive calls
   will return information on the next entry. The routine will return
   0 when there are no more entries.

   Header line: Must start with >. Will handle files which don't have
   the proper P1; or F1; parts of the header as well as those which
   do.

   Title line: Will read the name and source fields if correctly
   separated by a -, otherwise copies all information into the name.

   Sequence: May contain allowed puctuation. This will set the punct
   flag and information on the types found will be placed in seqinfo.
   White space and line breaks are ignored. Each chain should end with
   a *, but the routine will accept the last chain of an entry with no
   *. While the standard requires upper case text, this routine will
   handle lower case and convert it to upper case. While the routine
   does pretty well at last chains not terminated with a *, a last
   chain ending with a / not followed by a * but followed by a text
   line will be identified as incomplete rather than truncated.
   If the DoInsert flag is set, - signs in the sequence will be
   read as part of the sequence, otherwise they will be skipped. This
   is an addition to the PIR standard.

   Text lines: Text lines after an entry (beginning with R;, C;, A;, 
   N; or F;) are ignored.

   02.03.94 Original    By: ACRM
   03.03.94 Added / and = handling, upcasing, strcpy()->strncpy(),
            header lines without semi-colon, title lines without -
   07.03.94 Added sequence insertion handling and DoInsert parameter.
   11.05.94 buffer is now 504 characters (V38.0 spec allows 500 chars)
            Removes leading spaces from entry code and terminates at
            first space (V39.0 spec allows comments after the code).
   28.02.95 Added check that buffer doesn't overflow. Check on nseq
            changed to >=
   06.02.96 Removes trailing spaces from comment line
*/
int ReadPIR(FILE *fp, BOOL DoInsert, char **seqs, int maxchain, 
            SEQINFO *seqinfo, BOOL *punct, BOOL *error)
{
   int  ch,
        i,
        chpos,
        nseq = 0,
	ArraySize,
        SeqPos;
   char buffer[504],
        *ptr;
   BOOL InParen,
        GotStar;

   /* Initialise error and punct outputs                                */
   *error = FALSE;
   *punct = FALSE;

   /* Initialise seqinfo structure                                      */
   if(seqinfo != NULL)
   {
      seqinfo->code[0]    = '\0';
      seqinfo->name[0]    = '\0';
      seqinfo->source[0]  = '\0';
      seqinfo->fragment   = FALSE;
      seqinfo->paren      = FALSE;
      seqinfo->DotInParen = FALSE;
      seqinfo->NonExpJoin = FALSE;
      seqinfo->UnknownPos = FALSE;
      seqinfo->Incomplete = FALSE;
      seqinfo->Juxtapose  = FALSE;
      seqinfo->Truncated  = FALSE;
   }
   
   /* Skip over any characters until the first > sign                   */
   while((ch=fgetc(fp)) != EOF && ch != '>') ;

   /* Check for end of file                                             */
   if(ch==EOF) return(0);

   /* Read the rest of this line into a buffer                          */
   i = 0;
   while((ch=fgetc(fp)) != EOF && ch != '\n' && i<503)
      buffer[i++] = (char)ch;
   buffer[i] = '\0';

   /* Check for end of file                                             */
   if(ch==EOF) return(0);

   /* Set information in the seqinfo structure                          */
   if(seqinfo != NULL)
   {
      /* Fragment flag                                                  */
      if(buffer[2] == ';' && buffer[0] == 'F')
         seqinfo->fragment = TRUE;
      else 
         seqinfo->fragment = FALSE;

      /* Entry code                                                     */
      if(buffer[2] == ';')
      {
         KILLLEADSPACES(ptr,(buffer+3));
      }
      else
      {
         KILLLEADSPACES(ptr,buffer);
      }
      
      strncpy(seqinfo->code, ptr, 16);
      seqinfo->code[15] = '\0';

      /* Terminate entry code at first space since comments are allowed
         after the entry code (V39.0 spec)
      */
      for(i=0; seqinfo->code[i]; i++)
      {
         if(seqinfo->code[i] == ' ' || seqinfo->code[i] == '\t')
         {
            seqinfo->code[i] = '\0';
            break;
         }
      }
   }

   /* Now read the title line                                           */
   if(!fgets(buffer,240,fp))
      return(0);
   buffer[240] = '\0';

   /* 06.02.96 Remove any trailing spaces                               */
   KILLTRAILSPACES(buffer);

   /* Set information in the seqinfo structure                          */
   if(seqinfo)
   {
      TERMINATE(buffer);
      /* If it's a fully legal PIR file, there will be a - in the midle
         of the title line to separate name from source. If we don't
         find one, we copy the whole line into the name
      */
      if((ptr = strstr(buffer," - ")) != NULL)
      {
         *ptr = '\0';
         strncpy(seqinfo->source, ptr+3, 160);
         seqinfo->source[159] = '\0';
      }
      strncpy(seqinfo->name,   buffer, 160);
      seqinfo->name[159] = '\0';
      /* 06.02.96 Remove any trailing spaces                            */
      KILLTRAILSPACES(seqinfo->name);
   }

   /* Read the actual sequence info.                                    */
   chpos = 0;
   for(;;)
   {
      GotStar = FALSE;
      InParen = FALSE;

      /* Allocate some space for the sequence                           */
      ArraySize = ALLOCSIZE;
      if((seqs[nseq] = (char *)malloc(ArraySize * sizeof(char)))==NULL)
      {
         *error = TRUE;
         return(0);
      }

      SeqPos    = 0;
      
      /* Read characters, storing sequence and handling any 
         punctuation
      */
      while((ch = fgetc(fp)) != EOF && ch != '*' && ch != '>')
      {
         chpos++;
         
         if(isalpha(ch) || (ch == '-' && DoInsert))
         {
            /* This is a sequence entry (probably!)                     */
            seqs[nseq][SeqPos++] = (isupper(ch) ? ch : toupper(ch));

            /* If necessary, expand the sequence array                  */
            if(SeqPos >= ArraySize)
            {
               ArraySize += ALLOCSIZE;
               seqs[nseq] = (char *)realloc((void *)(seqs[nseq]), 
                                            ArraySize);
               if(seqs[nseq] == NULL)
               {
                  *error = TRUE;
                  return(0);
               }
            }
         }
         else if(ch == '/')
         {
            /* Sequence is incomplete or truncated                      */
            *punct = TRUE;
            
            if(seqinfo != NULL)
            {
               if(SeqPos == 0)   /* It's the first character in a chain */
               {

                  seqinfo->Truncated = TRUE;
               }
               else              /* Not first, is it last?              */
               {
                  /* Skip spaces and newlines till we get the next real
                     character
                  */
                  while((ch = fgetc(fp)) != EOF && 
                        (ch == ' ' || ch == '\t' || ch == '\n')) ;
                  /* Replace the character in the input stream          */
                  ungetc(ch,fp);
                  
                  if(ch == '*' || 
                     ch == EOF ||
                     ch == '>')               /* End of chain           */
                     seqinfo->Truncated  = TRUE;
                  else                        /* Middle of chain        */
                     seqinfo->Incomplete = TRUE;
               }
            }
         }
         else if(ch == '=')
         {
            /* Parts of the sequence may be juxtaposed                  */
            *punct = TRUE;
            if(seqinfo != NULL) seqinfo->Juxtapose = TRUE;
         }
         else if(ch == '(')
         {
            /* Start of a region in parentheses                         */
            InParen = TRUE;
            *punct = TRUE;
            if(seqinfo != NULL) seqinfo->paren = TRUE;
         }
         else if(ch == ')')
         {
            /* End of region in parentheses                             */
            InParen = FALSE;
            *punct = TRUE;
            if(seqinfo != NULL) seqinfo->paren = TRUE;
         }
         else if(ch == '.')
         {
            *punct = TRUE;

            if(InParen)
            {
               /* Previous aa >90% certain in position                  */
               if(seqinfo != NULL) seqinfo->DotInParen = TRUE;
            }
            else
            {
               /* Join in sequence not known experimentally but is clear
                  from sequence homology.
               */
               if(seqinfo != NULL) seqinfo->NonExpJoin = TRUE;
            }
         }
         else if(ch == ',')
         {
            /* Position of previous aa not known with confidence        */
            if(seqinfo != NULL) seqinfo->UnknownPos = TRUE;
         }
         else if(ch == '\n')
         {
            /* Start of new line, relevant to check on ;                */
            chpos = 0;
         }
         else if(ch == ';' && chpos == 2)
         {
            /* This is a text line, so the previous character wasn't
               a sequence item
            */
            SeqPos--;

            /* Ignore the rest of this line and reset chpos             */
            while((ch = fgetc(fp))!=EOF && ch != '\n') ;
            chpos = 0;
         }
      }  /* Reading this sequence                                       */

      /* Test the exit conditions from the read character loop          */
      if(ch == '*')
      {
         /* End of chain                                                */
         seqs[nseq][SeqPos] = '\0';
         GotStar = TRUE;
         if(++nseq >= maxchain)
         {
            *error = TRUE;
            return(nseq);
         }
      }
      else if(ch == '>')
      {
         /* Start of new entry                                          */
         ungetc(ch,fp);
         break;          /* Out of read for this sequence               */
      }
      else if(ch == EOF)
      {
         /* End of file                                                 */
         break;          /* Out of read for this sequence               */
      }
   }  /* Loop on with this sequence (next chain)                        */


   /* Now tidy up if we have an unfinished sequence                     */
   if(!GotStar) 
   {
      seqs[nseq][SeqPos] = '\0';
      if(!strlen(seqs[nseq]))
         free(seqs[nseq]);
      else
         nseq++;
   }

   return(nseq);
}


