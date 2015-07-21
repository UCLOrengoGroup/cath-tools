/*************************************************************************

   Program:    
   File:       align.c
   
   Version:    V3.3
   Date:       07.04.09
   Function:   Perform Needleman & Wunsch sequence alignment
   
   Copyright:  (c) SciTech Software 1993-2009
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
   A simple Needleman & Wunsch Dynamic Programming alignment of 2 
   sequences.  
   A window is not used so the routine may be a bit slow on long 
   sequences.

**************************************************************************

   Usage:
   ======
   First call ReadMDM() to read the mutation data matrix, then call
   align() to align the sequences.

**************************************************************************

   Revision History:
   =================
   V1.0  19.06.90 Original used in NW program
   V2.0  07.10.92 Original extracted from old NW program
   V2.1  16.06.93 Tidied for book
   V2.2  01.03.94 Changed static variable names
   V2.3  18.03.94 getc() -> fgetc()
   V2.4  24.11.94 ReadMDM() now looks after searching DATADIR
   V2.5  28.02.95 ReadMDM() and other code improved to cope with MDMs of
                  any size
   V2.6  26.07.95 Removed unused variables
   V2.7  21.08.95 Initialisation of matrix was incorrect leading to errors
                  at the end of the alignment.
   V2.8  24.08.95 calcscore() was doing an out-of-bounds array reference
                  if a character wasn't found
   V2.9  11.07.96 calcscore() changed to CalcMDMScore() and made 
                  non-static
   V2.10 09.09.96 Improved comments for ReadMDM()
   V2.11 17.09.96 Added ZeroMDM()
   V3.0  06.03.00 Traceback code rewritten to use a trace matrix created
                  while the main matrix is populated. New affinealign()
                  routine implemented. align() is now a wrapper to that.
   V3.1  06.02.03 Fixed for new version of GetWord()
   V3.2  27.02.07 Added affinealineuc() and CalcMDMScoreUC()
   V3.3  07.04.09 Complete re-write of ReadMDM() so it can read BLAST
                  style matrix files as well as our own

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "SysDefs.h"
#include "macros.h"
#include "array.h"
#include "general.h"
#include "seq.h"

/************************************************************************/
/* Defines and macros
*/
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define MAX3(c,d,e) (MAX(MAX((c),(d)),(e)))

#define DATAENV "DATADIR"   /* Environment variable or assign           */

#define MAXBUFF 400
#define MAXWORD 16

/* Type definition to store a X,Y coordinate pair in the matrix         */
typedef struct
{
   int x, y;
}  XY;


/************************************************************************/
/* Globals
*/
static int  **sMDMScore;
static char *sMDM_AAList;
static int  sMDMSize = 0;

/************************************************************************/
/* Prototypes
*/
static int  SearchForBest(int **matrix, int length1, int length2, 
                          int *BestI, int *BestJ, char *seq1, char *seq2, 
                          char *align1, char *align2);
static int  TraceBack(int **matrix, XY **dirn, int length1, int length2, 
                      char *seq1, char *seq2, char *align1, char *align2, 
                      int *align_len);


/************************************************************************/
/*>int align(char *seq1, int length1, char *seq2, int length2, 
             BOOL verbose, BOOL identity, int penalty,
             char *align1, char *align2, int *align_len)
   -----------------------------------------------------------
   Input:   char  *seq1         First sequence
            int   length1       First sequence length
            char  *seq2         Second sequence
            int   length2       Second sequence length
            BOOL  verbose       Display N&W matrix
            BOOL  identity      Use identity matrix
            int   penalty       Gap insertion penalty value
   Output:  char  *align1       Sequence 1 aligned
            char  *align2       Sequence 2 aligned
            int   *align_len    Alignment length
   Returns: int                 Alignment score (0 on error)
            
   Perform simple N&W alignment of seq1 and seq2. No window is used, so
   will be slow for long sequences.

   A single gap penalty is used, so gap extension incurrs no further
   penalty.

   Note that you must allocate sufficient memory for the aligned 
   sequences.
   The easy way to do this is to ensure that align1 and align2 are
   of length (length1+length2).

   06.03.00 Implemented as a wrapper to affinealign() which is the old
            align() routine, plus support for affine gap penalties,
            plus new traceback code based on storing the path as we
            go
*/
int align(char *seq1, 
          int  length1, 
          char *seq2, 
          int  length2, 
          BOOL verbose, 
          BOOL identity, 
          int  penalty, 
          char *align1, 
          char *align2,
          int  *align_len)
{
   return(affinealign(seq1, length1, seq2, length2, verbose, identity,
                      penalty, 0, align1, align2, align_len));
}


/************************************************************************/
/*>int affinealign(char *seq1, int length1, char *seq2, int length2, 
                   BOOL verbose, BOOL identity, int penalty, int penext, 
                   char *align1, char *align2, int *align_len)
   ---------------------------------------------------------------------
   Input:   char  *seq1         First sequence
            int   length1       First sequence length
            char  *seq2         Second sequence
            int   length2       Second sequence length
            BOOL  verbose       Display N&W matrix
            BOOL  identity      Use identity matrix
            int   penalty       Gap insertion penalty value
            int   penext        Extension penalty
   Output:  char  *align1       Sequence 1 aligned
            char  *align2       Sequence 2 aligned
            int   *align_len    Alignment length
   Returns: int                 Alignment score (0 on error)
            
   Perform simple N&W alignment of seq1 and seq2. No window is used, so
   will be slow for long sequences.

   Note that you must allocate sufficient memory for the aligned 
   sequences.
   The easy way to do this is to ensure that align1 and align2 are
   of length (length1+length2).

   07.10.92 Adapted from original written while at NIMR
   08.10.92 Split into separate routines
   09.10.92 Changed best structure to simple integers, moved 
            SearchForBest() into TraceBack()
   21.08.95 Was only filling in the bottom right cell at initialisation
            rather than all the right hand column and bottom row
   11.07.96 Changed calls to calcscore() to CalcMDMScore()
   06.03.00 Changed name to affinealign() (the routine align() is
            provided as a backwards compatible wrapper). Added penext 
            parameter. Now supports affine gap penalties with separate
            opening and extension penalties. The code now maintains
            the path as it goes.
**************************************************************************
******   NOTE AND CHANGES SHOULD BE PROPAGATED TO affinealignuc()   ******
**************************************************************************
*/
int affinealign(char *seq1, 
                int  length1, 
                char *seq2, 
                int  length2, 
                BOOL verbose, 
                BOOL identity, 
                int  penalty, 
                int  penext,
                char *align1, 
                char *align2,
                int  *align_len)
{
   XY    **dirn   = NULL;
   int   **matrix = NULL,
         maxdim,
         i,    j,    k,    l,
         i1,   j1,
         dia,  right, down,
         rcell, dcell, maxoff,
         match = 1,
         thisscore,
         gapext,
         score;
   
   maxdim = MAX(length1, length2);
   
   /* Initialise the score matrix                                       */
   if((matrix = (int **)Array2D(sizeof(int), maxdim, maxdim))==NULL)
      return(0);
   if((dirn   = (XY **)Array2D(sizeof(XY), maxdim, maxdim))==NULL)
      return(0);
      
   for(i=0;i<maxdim;i++)
   {
      for(j=0;j<maxdim;j++)
      {
         matrix[i][j] = 0;
         dirn[i][j].x = -1;
         dirn[i][j].y = -1;
      }
   }
    
   /* Fill in scores up the right hand side of the matrix               */
   for(j=0; j<length2; j++)
   {
      if(identity)
      {
         if(seq1[length1-1] == seq2[j]) matrix[length1-1][j] = match;
      }
      else
      {
         matrix[length1-1][j] = CalcMDMScore(seq1[length1-1], seq2[j]);
      }
   }

   /* Fill in scores along the bottom row of the matrix                 */
   for(i=0; i<length1; i++)
   {
      if(identity)
      {
         if(seq1[i] == seq2[length2-1]) matrix[i][length2-1] = match;
      }
      else
      {
         matrix[i][length2-1] = CalcMDMScore(seq1[i], seq2[length2-1]);
      }
   }

   i = length1 - 1;
   j = length2 - 1;
   
   /* Move back along the diagonal                                      */
   while(i > 0 && j > 0)
   {
      i--;
      j--;

      /* Fill in the scores along this row                              */
      for(i1 = i; i1 > -1; i1--)
      {
         dia   = matrix[i1+1][j+1];

         /* Find highest score to right of diagonal                     */
         rcell = i1+2;
         if(i1+2 >= length1)  right = 0;
         else                 right = matrix[i1+2][j+1] - penalty;
         
         gapext = 1;
         for(k = i1+3; k<length1; k++, gapext++)
         {
            thisscore = matrix[k][j+1] - (penalty + gapext*penext);
            
            if(thisscore > right) 
            {
               right = thisscore;
               rcell = k;
            }
         }

         /* Find highest score below diagonal                           */
         dcell = j+2;
         if(j+2 >= length2)  down = 0;
         else                down   = matrix[i1+1][j+2] - penalty;
         
         gapext = 1;
         for(l = j+3; l<length2; l++, gapext++)
         {
            thisscore = matrix[i1+1][l] - (penalty + gapext*penext);

            if(thisscore > down) 
            {
               down = thisscore;
               dcell = l;
            }
         }
         
         /* Set score to best of these                                  */
         maxoff = MAX(right, down);
         if(dia >= maxoff)
         {
            matrix[i1][j] = dia;
            dirn[i1][j].x = i1+1;
            dirn[i1][j].y = j+1;
         }
         else
         {
            if(right > down)
            {
               matrix[i1][j] = right;
               dirn[i1][j].x = rcell;
               dirn[i1][j].y = j+1;
            }
            else
            {
               matrix[i1][j] = down;
               dirn[i1][j].x = i1+1;
               dirn[i1][j].y = dcell;
            }
         }
       
         /* Add the score for a match                                   */
         if(identity)
         {
            if(seq1[i1] == seq2[j]) matrix[i1][j] += match;
         }
         else
         {
            matrix[i1][j] += CalcMDMScore(seq1[i1],seq2[j]);
         }
      }

      /* Fill in the scores in this column                              */
      for(j1 = j; j1 > -1; j1--)
      {
         dia   = matrix[i+1][j1+1];
         
         /* Find highest score to right of diagonal                     */
         rcell = i+2;
         if(i+2 >= length1)   right = 0;
         else                 right = matrix[i+2][j1+1] - penalty;

         gapext = 1;
         for(k = i+3; k<length1; k++, gapext++)
         {
            thisscore = matrix[k][j1+1] - (penalty + gapext*penext);
            
            if(thisscore > right) 
            {
               right = thisscore;
               rcell = k;
            }
         }

         /* Find highest score below diagonal                           */
         dcell = j1+2;
         if(j1+2 >= length2)  down = 0;
         else                 down = matrix[i+1][j1+2] - penalty;

         gapext = 1;
         for(l = j1+3; l<length2; l++, gapext++)
         {
            thisscore = matrix[i+1][l] - (penalty + gapext*penext);
            
            if(thisscore > down) 
            {
               down = thisscore;
               dcell = l;
            }
         }

         /* Set score to best of these                                  */
         maxoff = MAX(right, down);
         if(dia >= maxoff)
         {
            matrix[i][j1] = dia;
            dirn[i][j1].x = i+1;
            dirn[i][j1].y = j1+1;
         }
         else
         {
            if(right > down)
            {
               matrix[i][j1] = right;
               dirn[i][j1].x = rcell;
               dirn[i][j1].y = j1+1;
            }
            else
            {
               matrix[i][j1] = down;
               dirn[i][j1].x = i+1;
               dirn[i][j1].y = dcell;
            }
         }
       
         /* Add the score for a match                                   */
         if(identity)
         {
            if(seq1[i] == seq2[j1]) matrix[i][j1] += match;
         }
         else
         {
            matrix[i][j1] += CalcMDMScore(seq1[i],seq2[j1]);
         }
      }
   } 
   
   score = TraceBack(matrix, dirn, length1, length2,
                     seq1, seq2, align1, align2, align_len);

   if(verbose)
   {
      printf("Matrix:\n-------\n");
      for(j=0; j<length2;j++)
      {
         for(i=0; i<length1; i++)
         {
            printf("%3d ",matrix[i][j]);
         }
         printf("\n");
      }

      printf("Path:\n-----\n");
      for(j=0; j<length2;j++)
      {
         for(i=0; i<length1; i++)
         {
            printf("(%3d,%3d) ",dirn[i][j].x,dirn[i][j].y);
         }
         printf("\n");
      }
   }
    
   FreeArray2D((char **)matrix, maxdim, maxdim);
   FreeArray2D((char **)dirn,   maxdim, maxdim);
    
   return(score);
}


/************************************************************************/
/*>int affinealign(char *seq1, int length1, char *seq2, int length2, 
                   BOOL verbose, BOOL identity, int penalty, int penext, 
                   char *align1, char *align2, int *align_len)
   ---------------------------------------------------------------------
   Input:   char  *seq1         First sequence
            int   length1       First sequence length
            char  *seq2         Second sequence
            int   length2       Second sequence length
            BOOL  verbose       Display N&W matrix
            BOOL  identity      Use identity matrix
            int   penalty       Gap insertion penalty value
            int   penext        Extension penalty
   Output:  char  *align1       Sequence 1 aligned
            char  *align2       Sequence 2 aligned
            int   *align_len    Alignment length
   Returns: int                 Alignment score (0 on error)
            
   Perform simple N&W alignment of seq1 and seq2. No window is used, so
   will be slow for long sequences.

   Note that you must allocate sufficient memory for the aligned 
   sequences.
   The easy way to do this is to ensure that align1 and align2 are
   of length (length1+length2).

   07.10.92 Adapted from original written while at NIMR
   08.10.92 Split into separate routines
   09.10.92 Changed best structure to simple integers, moved 
            SearchForBest() into TraceBack()
   21.08.95 Was only filling in the bottom right cell at initialisation
            rather than all the right hand column and bottom row
   11.07.96 Changed calls to calcscore() to CalcMDMScore()
   06.03.00 Changed name to affinealign() (the routine align() is
            provided as a backwards compatible wrapper). Added penext 
            parameter. Now supports affine gap penalties with separate
            opening and extension penalties. The code now maintains
            the path as it goes.
   27.02.07 Exactly as affinealign() but upcases characters before
            comparison

**************************************************************************
******    NOTE AND CHANGES SHOULD BE PROPAGATED TO affinealign()    ******
**************************************************************************
*/
int affinealignuc(char *seq1, 
                  int  length1, 
                  char *seq2, 
                  int  length2, 
                  BOOL verbose, 
                  BOOL identity, 
                  int  penalty, 
                  int  penext,
                  char *align1, 
                  char *align2,
                  int  *align_len)
{
   XY    **dirn   = NULL;
   int   **matrix = NULL,
         maxdim,
         i,    j,    k,    l,
         i1,   j1,
         dia,  right, down,
         rcell, dcell, maxoff,
         match = 1,
         thisscore,
         gapext,
         score;
   
   maxdim = MAX(length1, length2);
   
   /* Initialise the score matrix                                       */
   if((matrix = (int **)Array2D(sizeof(int), maxdim, maxdim))==NULL)
      return(0);
   if((dirn   = (XY **)Array2D(sizeof(XY), maxdim, maxdim))==NULL)
      return(0);
      
   for(i=0;i<maxdim;i++)
   {
      for(j=0;j<maxdim;j++)
      {
         matrix[i][j] = 0;
         dirn[i][j].x = -1;
         dirn[i][j].y = -1;
      }
   }
    
   /* Fill in scores up the right hand side of the matrix               */
   for(j=0; j<length2; j++)
   {
      if(identity)
      {
         if(seq1[length1-1] == seq2[j]) matrix[length1-1][j] = match;
      }
      else
      {
         matrix[length1-1][j] = CalcMDMScoreUC(seq1[length1-1], seq2[j]);
      }
   }

   /* Fill in scores along the bottom row of the matrix                 */
   for(i=0; i<length1; i++)
   {
      if(identity)
      {
         if(seq1[i] == seq2[length2-1]) matrix[i][length2-1] = match;
      }
      else
      {
         matrix[i][length2-1] = CalcMDMScoreUC(seq1[i], seq2[length2-1]);
      }
   }

   i = length1 - 1;
   j = length2 - 1;
   
   /* Move back along the diagonal                                      */
   while(i > 0 && j > 0)
   {
      i--;
      j--;

      /* Fill in the scores along this row                              */
      for(i1 = i; i1 > -1; i1--)
      {
         dia   = matrix[i1+1][j+1];

         /* Find highest score to right of diagonal                     */
         rcell = i1+2;
         if(i1+2 >= length1)  right = 0;
         else                 right = matrix[i1+2][j+1] - penalty;
         
         gapext = 1;
         for(k = i1+3; k<length1; k++, gapext++)
         {
            thisscore = matrix[k][j+1] - (penalty + gapext*penext);
            
            if(thisscore > right) 
            {
               right = thisscore;
               rcell = k;
            }
         }

         /* Find highest score below diagonal                           */
         dcell = j+2;
         if(j+2 >= length2)  down = 0;
         else                down   = matrix[i1+1][j+2] - penalty;
         
         gapext = 1;
         for(l = j+3; l<length2; l++, gapext++)
         {
            thisscore = matrix[i1+1][l] - (penalty + gapext*penext);

            if(thisscore > down) 
            {
               down = thisscore;
               dcell = l;
            }
         }
         
         /* Set score to best of these                                  */
         maxoff = MAX(right, down);
         if(dia >= maxoff)
         {
            matrix[i1][j] = dia;
            dirn[i1][j].x = i1+1;
            dirn[i1][j].y = j+1;
         }
         else
         {
            if(right > down)
            {
               matrix[i1][j] = right;
               dirn[i1][j].x = rcell;
               dirn[i1][j].y = j+1;
            }
            else
            {
               matrix[i1][j] = down;
               dirn[i1][j].x = i1+1;
               dirn[i1][j].y = dcell;
            }
         }
       
         /* Add the score for a match                                   */
         if(identity)
         {
            if(seq1[i1] == seq2[j]) matrix[i1][j] += match;
         }
         else
         {
            matrix[i1][j] += CalcMDMScoreUC(seq1[i1],seq2[j]);
         }
      }

      /* Fill in the scores in this column                              */
      for(j1 = j; j1 > -1; j1--)
      {
         dia   = matrix[i+1][j1+1];
         
         /* Find highest score to right of diagonal                     */
         rcell = i+2;
         if(i+2 >= length1)   right = 0;
         else                 right = matrix[i+2][j1+1] - penalty;

         gapext = 1;
         for(k = i+3; k<length1; k++, gapext++)
         {
            thisscore = matrix[k][j1+1] - (penalty + gapext*penext);
            
            if(thisscore > right) 
            {
               right = thisscore;
               rcell = k;
            }
         }

         /* Find highest score below diagonal                           */
         dcell = j1+2;
         if(j1+2 >= length2)  down = 0;
         else                 down = matrix[i+1][j1+2] - penalty;

         gapext = 1;
         for(l = j1+3; l<length2; l++, gapext++)
         {
            thisscore = matrix[i+1][l] - (penalty + gapext*penext);
            
            if(thisscore > down) 
            {
               down = thisscore;
               dcell = l;
            }
         }

         /* Set score to best of these                                  */
         maxoff = MAX(right, down);
         if(dia >= maxoff)
         {
            matrix[i][j1] = dia;
            dirn[i][j1].x = i+1;
            dirn[i][j1].y = j1+1;
         }
         else
         {
            if(right > down)
            {
               matrix[i][j1] = right;
               dirn[i][j1].x = rcell;
               dirn[i][j1].y = j1+1;
            }
            else
            {
               matrix[i][j1] = down;
               dirn[i][j1].x = i+1;
               dirn[i][j1].y = dcell;
            }
         }
       
         /* Add the score for a match                                   */
         if(identity)
         {
            if(seq1[i] == seq2[j1]) matrix[i][j1] += match;
         }
         else
         {
            matrix[i][j1] += CalcMDMScoreUC(seq1[i],seq2[j1]);
         }
      }
   } 
   
   score = TraceBack(matrix, dirn, length1, length2,
                     seq1, seq2, align1, align2, align_len);

   if(verbose)
   {
      printf("Matrix:\n-------\n");
      for(j=0; j<length2;j++)
      {
         for(i=0; i<length1; i++)
         {
            printf("%3d ",matrix[i][j]);
         }
         printf("\n");
      }

      printf("Path:\n-----\n");
      for(j=0; j<length2;j++)
      {
         for(i=0; i<length1; i++)
         {
            printf("(%3d,%3d) ",dirn[i][j].x,dirn[i][j].y);
         }
         printf("\n");
      }
   }
    
   FreeArray2D((char **)matrix, maxdim, maxdim);
   FreeArray2D((char **)dirn,   maxdim, maxdim);
    
   return(score);
}


/************************************************************************/
/*>BOOL ReadMDM(char *mdmfile)
   ---------------------------
   Input:   char  *mdmfile    Mutation data matrix filename
   Returns: BOOL              Success?
   
   Read mutation data matrix into static global arrays. The matrix may
   have comments at the start introduced with a ! in the first column.
   The matrix must be complete (i.e. a triangular matrix will not
   work). A line describing the residue types must appear, and may
   be placed before or after the matrix itself

   07.10.92 Original
   18.03.94 getc() -> fgetc()
   24.11.94 Automatically looks in DATAENV if not found in current 
            directory
   28.02.95 Modified to read any size MDM and allow comments
            Also allows the list of aa types before or after the actual
            matrix
   26.07.95 Removed unused variables
   06.02.03 Fixed for new version of GetWord()
   07.04.09 Completely re-written to allow it to read BLAST style matrix
            files as well as the ones used previously
            Allow comments introduced with # as well as !
            Uses MAXWORD rather than hardcoded 16
*/
BOOL ReadMDM(char *mdmfile)
{
   FILE *mdm = NULL;
   int  i, j, k, row, tmpStoreSize;
   char buffer[MAXBUFF],
        word[MAXWORD],
        *p,
        **tmpStore;
   BOOL noenv;

   if((mdm=OpenFile(mdmfile, DATAENV, "r", &noenv))==NULL)
   {
      return(FALSE);
   }

   /* First read the file to determine the dimensions                   */
   while(fgets(buffer,MAXBUFF,mdm))
   {
      TERMINATE(buffer);
      KILLLEADSPACES(p,buffer);

      /* First line which is non-blank and non-comment                  */
      if(strlen(p) && p[0] != '!' && p[0] != '#')
      {
         sMDMSize = 0;
         for(p = buffer; p!=NULL;)
         {
            p = GetWord(p, word, MAXWORD);
            /* Increment counter if this is numeric                     */
            if(isdigit(word[0]) || 
               ((word[0] == '-')&&(isdigit(word[1]))))
               sMDMSize++;
         }
         if(sMDMSize)
            break;
      }
   }

   /* Allocate memory for the MDM and the AA List                       */
   if((sMDMScore = (int **)Array2D(sizeof(int),sMDMSize,sMDMSize))==NULL)
      return(FALSE);
   if((sMDM_AAList = (char *)malloc((sMDMSize+1)*sizeof(char)))==NULL)
   {
      FreeArray2D((char **)sMDMScore, sMDMSize, sMDMSize);
      return(FALSE);
   }

   /* Allocate temporary storage for a row from the matrix              */
   tmpStoreSize = 2*sMDMSize;
   if((tmpStore = (char **)Array2D(sizeof(char), tmpStoreSize, MAXWORD))
      ==NULL)
   {
      free(sMDM_AAList);
      FreeArray2D((char **)sMDMScore, sMDMSize, sMDMSize);
      return(FALSE);
   }

   /* Fill the matrix with zeros                                        */
   for(i=0; i<sMDMSize; i++)
   {
      for(j=0; j<sMDMSize; j++)
      {
         sMDMScore[i][j] = 0;
      }
   }

   /* Rewind the file and read the actual data                          */
   rewind(mdm);
   row = 0;
   while(fgets(buffer,MAXBUFF,mdm))
   {
      int Numeric;
      
      TERMINATE(buffer);
      KILLLEADSPACES(p,buffer);

      /* Check line is non-blank and non-comment                        */
      if(strlen(p) && p[0] != '!' && p[0] != '#')
      {
         Numeric = 0;
         for(p = buffer, i = 0; p!=NULL && i<tmpStoreSize; i++)
         {
            p = GetWord(p, tmpStore[i], MAXWORD);
            /* Incremement Numeric counter if it's a numeric field      */
            if(isdigit(tmpStore[i][0]) || 
               ((tmpStore[i][0] == '-')&&(isdigit(tmpStore[i][1]))))
            {
               Numeric++;
            }
         }

         /* No numeric fields so it is the amino acid names             */
         if(Numeric == 0)
         {
            for(j = 0; j<i && j<sMDMSize; j++)
            {
               sMDM_AAList[j] = tmpStore[j][0];
            }
         }
         else
         {
            /* There were numeric fields, so copy them into the matrix,
               skipping any non-numeric fields
               j counts the input fields
               k counts the fields in sMDMScore
               row counts the row in sMDMScore
            */
            for(j=0, k=0; j<i && k<sMDMSize; j++)
            {
               if(isdigit(tmpStore[j][0]) || 
                  ((tmpStore[j][0] == '-')&&(isdigit(tmpStore[j][1]))))
               {
                  sscanf(tmpStore[j],"%d",&(sMDMScore[row][k]));
                  k++;
               }
            }
            
            row++;
         }
      }
   }
   fclose(mdm);
   FreeArray2D((char **)tmpStore, tmpStoreSize, MAXWORD);
   
   return(TRUE);
}


/************************************************************************/
/*>static int SearchForBest(int **matrix, int length1, int length2, 
                            int *BestI, int *BestJ, char *seq1, 
                            char *seq2, char *align1, char *align2)
   ----------------------------------------------------------------
   Input:   int  **matrix   N&W matrix
            int  length1    Length of first sequence
            int  length2    Length of second sequence
            int  *BestI     x position of highest score
            int  *BestJ     y position of highest score
            char *seq1      First sequence
            char *seq2      Second sequence
   Output:  char *align1    First sequence with end aligned correctly
            char *align2    Second sequence with end aligned correctly
   Returns: int             Alignment length thus far

   Searches the outside of the matrix for the best score and starts the
   alignment by putting in any starting - characters.

   08.10.92 Original extracted from Align()
*/
static int SearchForBest(int  **matrix, 
                         int  length1, 
                         int  length2, 
                         int  *BestI, 
                         int  *BestJ,
                         char *seq1, 
                         char *seq2, 
                         char *align1, 
                         char *align2)
{
   int   ai, 
         besti,   bestj, 
         i,       j;
   
   /* Now search the outside of the matrix for the highest scoring cell */
   ai    = 0;
   besti = 0;
   for(i = 1; i < length1; i++) 
   {
      if(matrix[i][0] > matrix[besti][0]) besti = i;
   }
   bestj = 0;
   for(j = 1; j < length2; j++)
   {
      if(matrix[0][j] > matrix[0][bestj]) bestj = j;
   }
   if(matrix[besti][0] > matrix[0][bestj])
   {
      *BestI = besti;
      *BestJ = 0;
      for(i=0; i<*BestI; i++)
      {
         align1[ai] = seq1[i];
         align2[ai++] = '-';
      }
   }
   else
   {
      *BestI = 0;
      *BestJ = bestj;
      for(j=0; j<*BestJ; j++)
      {
         align1[ai] = '-';
         align2[ai++] = seq2[j];
      }
   }
   return(ai);
}


/************************************************************************/
/*>static int TraceBack(int **matrix, XY **dirn, 
                        int length1, int length2, 
                        char *seq1, char *seq2, char *align1, 
                        char *align2, int *align_len)
   ----------------------------------------------------------
   Input:   int  **matrix   N&W matrix
            XY   **dirn     Direction Matrix
            int  length1    Length of first sequence
            int  length2    Length of second sequence
            char *seq1      First sequence
            char *seq2      Second sequence
   Output:  char *align1    First sequence aligned
            char *align2    Second sequence aligned
            int  *align_len Aligned sequence length
   Returns: int             Alignment score

   Does the traceback to find the aligment.

   08.10.92 Original extracted from Align(). Rewritten to do tracing
            correctly.
   09.10.92 Changed to call SearchForBest(). Nor returns score rather than
            length.
   06.03.00 Recoded to take the path matrix which is calculated within
            the main affinealign() routine rather than calculating the
            path as we go. penalty parameter removed as this is no longer
            needed. dirn parameter added.
   28.09.00 Fixed bug where last inserts were printed properly if one chain
            ended first
*/
static int TraceBack(int  **matrix, 
                     XY   **dirn,
                     int  length1, 
                     int  length2, 
                     char *seq1, 
                     char *seq2, 
                     char *align1, 
                     char *align2, 
                     int  *align_len)
{
   int   i,    j, 
         ai, 
         BestI,BestJ;
   XY    nextCell;

   ai = SearchForBest(matrix, length1, length2, &BestI, &BestJ, 
                       seq1, seq2, align1, align2);

   /* Now trace back to find the alignment                              */
   i            = BestI;
   j            = BestJ;
   align1[ai]   = seq1[i];
   align2[ai++] = seq2[j];

   while(i < length1-1 && j < length2-1)
   {
      nextCell.x = dirn[i][j].x;
      nextCell.y = dirn[i][j].y;
      if((nextCell.x == i+1) && (nextCell.y == j+1))
      {
         /* We are inheriting from the diagonal                         */
         i++;
         j++;
      }
      else if(nextCell.y == j+1)
      {
         /* We are inheriting from the off-diagonal inserting a gap in
            the y-sequence (seq2)
         */
         i++;
         j++;
         while((i < nextCell.x) && (i < length1-1))
         {
            align1[ai] = seq1[i++];
            align2[ai++] = '-';
         }
      }
      else if(nextCell.x == i+1)
      {
         /* We are inheriting from the off-diagonal inserting a gap in
            the x-sequence (seq1)
         */
         i++;
         j++;
         while((j < nextCell.y) && (j < length2-1))
         {
            align1[ai] = '-';
            align2[ai++] = seq2[j++];
         }
      }
      else
      {
         /* Cockup!                                                     */
         fprintf(stderr,"align.c/TraceBack() internal error\n");
      }
      
      align1[ai]   = seq1[i];
      align2[ai++] = seq2[j];
   }

   /* If one sequence finished first, fill in the end with insertions   */
   if(i < length1-1)
   {
      for(j=i+1; j<length1; j++)
      {
         align1[ai]   = seq1[j];
         align2[ai++] = '-';
      }
   }
   else if(j < length2-1)
   {
      for(i=j+1; i<length2; i++)
      {
         align1[ai]   = '-';
         align2[ai++] = seq2[i];
      }
   }
   
   *align_len = ai;
   
   return(matrix[BestI][BestJ]);
}


/************************************************************************/
/*>int CalcMDMScore(char resa, char resb)
   --------------------------------------
   Input:   char   resa      First residue
            char   resb      Second residue
   Returns: int              score

   Calculate score from static globally stored mutation data matrix

   07.10.92 Adapted from NIMR-written original
   24.11.94 Only gives 10 warnings
   28.02.95 Modified to use sMDMSize
   24.08.95 If a residue was not found was doing an out-of-bounds array
            reference causing a potential core dump
   11.07.96 Name changed from calcscore() and now non-static
*/
int CalcMDMScore(char resa, char resb)
{
   int        i,j;
   static int NWarn = 0;
   BOOL       Warned = FALSE;

   for(i=0;i<sMDMSize;i++)
   {
      if(resa==sMDM_AAList[i]) break;
   }
   if(i==sMDMSize) 
   {
      if(NWarn < 10)
         printf("Residue %c not found in matrix\n",resa);
      else if(NWarn == 10)
         printf("More residues not found in matrix...\n");
      Warned = TRUE;
   }
   for(j=0;j<sMDMSize;j++)
   {
      if(resb==sMDM_AAList[j]) break;
   }
   if(j==sMDMSize) 
   {
      if(NWarn < 10)
         printf("Residue %c not found in matrix\n",resb);
      else if(NWarn == 10)
         printf("More residues not found in matrix...\n");
      Warned = TRUE;
   }
   
   if(Warned)
   { 
      NWarn++;
      return(0);
   }

   return(sMDMScore[i][j]);
}                               

/************************************************************************/
/*>int CalcMDMScoreUC(char resa, char resb)
   ----------------------------------------
   Input:   char   resa      First residue
            char   resb      Second residue
   Returns: int              score

   Calculate score from static globally stored mutation data matrix

   07.10.92 Adapted from NIMR-written original
   24.11.94 Only gives 10 warnings
   28.02.95 Modified to use sMDMSize
   24.08.95 If a residue was not found was doing an out-of-bounds array
            reference causing a potential core dump
   11.07.96 Name changed from calcscore() and now non-static
   27.02.07 As CalcMDMScore() but upcases characters before comparison
*/
int CalcMDMScoreUC(char resa, char resb)
{
   int        i,j;
   static int NWarn = 0;
   BOOL       Warned = FALSE;

   resa = (islower(resa)?toupper(resa):resa);
   resb = (islower(resb)?toupper(resb):resb);

   for(i=0;i<sMDMSize;i++)
   {
      if(resa==sMDM_AAList[i]) break;
   }
   if(i==sMDMSize) 
   {
      if(NWarn < 10)
         printf("Residue %c not found in matrix\n",resa);
      else if(NWarn == 10)
         printf("More residues not found in matrix...\n");
      Warned = TRUE;
   }
   for(j=0;j<sMDMSize;j++)
   {
      if(resb==sMDM_AAList[j]) break;
   }
   if(j==sMDMSize) 
   {
      if(NWarn < 10)
         printf("Residue %c not found in matrix\n",resb);
      else if(NWarn == 10)
         printf("More residues not found in matrix...\n");
      Warned = TRUE;
   }
   
   if(Warned)
   { 
      NWarn++;
      return(0);
   }

   return(sMDMScore[i][j]);
}                               

/************************************************************************/
/*>int ZeroMDM(void)
   -----------------
   Returns: int           Maximum value in modified matrix

   Modifies all values in the MDM such that the minimum value is 0
   17.09.96 Original
*/
int ZeroMDM(void)
{
   int MinVal = sMDMScore[0][0],
       MaxVal = sMDMScore[0][0],
       i, j;

   /* Find the minimum and maximum values on the matrix                 */
   for(i=0; i<sMDMSize; i++)
   {
      for(j=0; j<sMDMSize; j++)
      {
         if(sMDMScore[i][j] < MinVal)
         {
            MinVal = sMDMScore[i][j];
         }
         else if(sMDMScore[i][j] > MaxVal)
         {
            MaxVal = sMDMScore[i][j];
         }
      }
   }
   
   /* Now subtract the MinVal from all cells in the matrix so it starts
      at zero.
   */
   for(i=0; i<sMDMSize; i++)
   {
      for(j=0; j<sMDMSize; j++)
      {
         sMDMScore[i][j] -= MinVal;
      }
   }
   
   /* Return maximum value in modified matrix                           */
   return(MaxVal-MinVal);
}

               
            
      
#ifdef DEMO   
int main(int argc, char **argv)
{
   char seq1[] = "ACTCLMCT",
        seq2[] = "ACTCCT",
        align1[100],
        align2[100];
   int  score, al_len;
   int  i, j;
   
   ReadMDM("pet91.mat");

   for(i=0; i<sMDMSize; i++)
   {
      printf("  %c", sMDM_AAList[i]);
   }
   printf("\n");
   
   for(i=0; i<sMDMSize; i++)
   {
      for(j=0; j<sMDMSize; j++)
      {
         printf("%3d", sMDMScore[i][j]);
      }
      printf("\n");
   }
   
   score = affinealign(seq1, strlen(seq1), seq2, strlen(seq2), 
                       TRUE, FALSE,
                       10, 1, align1, align2, &al_len);

   align1[al_len] = '\0';
   align2[al_len] = '\0';

   printf("%s\n", align1);
   printf("%s\n", align2);
   
   return(0);
}
#endif
