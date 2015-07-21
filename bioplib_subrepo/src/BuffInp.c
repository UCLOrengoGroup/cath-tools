/*************************************************************************

   Program:    
   File:       BuffInp.c
   
   Version:    V1.1R
   Date:       11.03.94
   Function:   Read from a file a line at a time, allowing one to probe
               ahead and look at the contants of the next line without
               removing it from the input stream.
   
   Copyright:  (c) SciTech Software 1994
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
   V1.0  08.03.94 Original
   V1.1  11.03.94 Changes to ReadBufferedFile() and ProbeBufferedFile()
                  to RePrompt() if reading from stdin when we get a 
                  blank line

*************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include <string.h>

#include "macros.h"
#include "BuffInp.h"
#include "WindIO.h"

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
/*>INBUFFER OpenBufferedFile(char *filename, int maxstr)
   -----------------------------------------------------
   Input:   char     *filename   File name
            int      maxstr      Max string length in file for buffering
   Returns: INBUFFER *           Pointer to a buffered file

   Open a file for buffered input. This allows probe-ahead to look at the
   contents of the next line without removing it from the input stream.

   28.02.94 Original    By: ACRM
   03.03.94 If filename is NULL, make file stdin
*/
INBUFFER *OpenBufferedFile(char *filename, int maxstr)
{
   FILE     *fp;
   INBUFFER *BuffStruc = NULL;
   
   if(filename != NULL)
      fp=fopen(filename,"r");
   else
      fp=stdin;

   if(fp!=NULL)
   {
      if((BuffStruc = (INBUFFER *)malloc(sizeof(INBUFFER)))!=NULL)
      {
         BuffStruc->fp         = fp;
         BuffStruc->nlines     = 0;
         BuffStruc->maxstr     = maxstr;
         if((BuffStruc->buffer = 
             (char *)malloc(maxstr * sizeof(char))) != NULL)
            return(BuffStruc);
      }
      fclose(fp);
   }
   return(NULL);
}

/************************************************************************/
/*>BOOL ReadBufferedFile(INBUFFER *bfp, char *string, int length)
   --------------------------------------------------------------
   Input:   INBUFFER *bfp      Pointer to a buffered file structure
            int      length    Size of output string
   Output:  char     *string   Output string read from file
   Returns: BOOL               TRUE: Successful read
                               FALSE: End of file (or error)

   Reads a line from a buffered file (like fgets()).
   Blank lines in the file will be skipped.

   28.02.94 Original    By: ACRM
   07.03.94 Added code to skip blank lines
   10.03.94 Added call to RePrompt() if we're reading from stdin and
            we get a blank line
*/
BOOL ReadBufferedFile(INBUFFER *bfp, char *string, int length)
{
   int bufflen = 0;
   
   if(bfp == NULL)
      return(FALSE);

   if(bfp->nlines == 0)
   {
      while(bufflen==0)
      {
         if(fgets(string,length,bfp->fp))
         {
            TERMINATE(string);
            bufflen=strlen(string);
         }
         else
         {
            return(FALSE);
         }

         if(!bufflen && bfp->fp == stdin)
            RePrompt();
      }
   }
   else
   {
      strncpy(string,bfp->buffer,length);
      (bfp->nlines)--;
   }

   return(TRUE);
}

/************************************************************************/
/*>BOOL ProbeBufferedFile(INBUFFER *bfp, char *string, int length)
   ---------------------------------------------------------------
   Input:   INBUFFER *bfp      Pointer to a buffered file structure
            int      length    Size of output string
   Output:  char     *string   Output string read from file
   Returns: BOOL               TRUE: Successful read
                               FALSE: End of file (or error)

   Read the next line from a buffered file without removing it from
   the input stream. Repeated calls will thus return the same string.
   The next call to ReadBufferedFile will also output the same string,
   but will remove the line from the input stream.
   Blank lines in the file will be skipped.

   28.02.94 Original    By: ACRM
   07.03.94 Added code to skip blank lines
   10.03.94 Added call to RePrompt() after blank line
*/
BOOL ProbeBufferedFile(INBUFFER *bfp, char *string, int length)
{
   int bufflen = 0;
   
   if(bfp==NULL)
      return(FALSE);

   while(bufflen==0)
   {
      if(bfp->nlines == 0)
      {
         if(fgets(bfp->buffer,bfp->maxstr,bfp->fp))
         {
            TERMINATE(bfp->buffer);

            bufflen = strlen(bfp->buffer);
            if(bufflen)   (bfp->nlines)++;
         }
         else
         {
            return(FALSE);
         }

         if(!bufflen && bfp->fp == stdin)
            RePrompt();
      }
   }

   strncpy(string,bfp->buffer,length);
   
   return(TRUE);
}

