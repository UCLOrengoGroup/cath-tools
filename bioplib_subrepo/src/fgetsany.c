/*************************************************************************

   Program:    
   File:       fgetsany.c
   
   Version:    V1.0
   Date:       11.09.95
   Function:   Like fgets(), but allocates memory and returns pointer
               to memory block
   
   Copyright:  (c) SciTech Software 1995
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

   fgetsany() provides a routine like fgets() for reading strings from
   a file, but does not require you to impose a limit on the length of
   string which may be read. fgetsany() allocates memory to store a
   string of any length and returns a pointer to that allocated memory.
   After use, this memory must be freed by the calling routine.

**************************************************************************

   Usage:
   ======

   fgetsany() returns NULL on end-of-file and if memory allocation
   failed. It uses the global `errno' variable to indicate a memory
   allocation failure (errno==ENOMEM). 

   Typical usage is as follows:



   #include "bioplib/general.h"

   main()
   {
      char *ptr;
      FILE *fp;
   
      fp = fopen("test.txt","r");
   
      while((ptr=fgetsany(fp))!=NULL)
      {
         printf("%s",ptr);
         free(ptr);
      }
   
      if(errno)
         perror("SYSTEM ERROR");
   }


**************************************************************************

   Revision History:
   =================

   V1.0  11.09.95 Original    By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 80   /* Memory allocated in chunks of this size         */

/************************************************************************/
/* Globals
*/

extern int errno;

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>char *fgetsany(FILE *fp)
   ------------------------
   Input:   FILE    *fp           File pointer open for reading
   Returns: char    *             Allocated string or NULL on EOF and
                                  no memory (errno==ENOMEM)

   fgetsany() provides a routine like fgets() for reading strings from
   a file, but does not require you to impose a limit on the length of
   string which may be read. fgetsany() allocates memory to store a
   string of any length and returns a pointer to that allocated memory.
   After use, this memory must be freed by the calling routine.

   fgetsany() returns NULL on end-of-file and if memory allocation
   failed. It uses the global `errno' variable to indicate a memory
   allocation failure (errno==ENOMEM). 

   11.09.95 Original    By: ACRM
*/
char *fgetsany(FILE *fp)
{
   int  bufflen = MAXBUFF,
        buffpos = 0,
        ch;
   char *buffer = NULL;

   /* Get a character and test for end of file                          */
   if((ch = fgetc(fp))==EOF)
      return(NULL);

   /* If it wasn't EOF, put the character back in the input stream      */
   ungetc(ch, fp);
   
   /* OK, wasn't EOF, so allocate initial buffer space                  */
   if((buffer = (char *)malloc(bufflen * sizeof(char)))==NULL)
   {
      errno = ENOMEM;
      return(NULL);
   }

   /* Read characters from file                                         */
   do 
   {
      ch = fgetc(fp);

      /* If it's not end of file, store the character                   */
      if(ch != EOF)
      {
         buffer[buffpos] = ch;

         /* Increment buffer position and test for overflow             */
         if(++buffpos == bufflen)
         {
            bufflen += MAXBUFF;
   
            /* Allocate additional space on overflow                    */
            if((buffer = realloc(buffer, bufflen))==NULL)
            {
               errno = ENOMEM;
               return(NULL);
            }
         }
      }
   }  while((ch != EOF) && (ch != '\n'));
   
   /* Terminate the string                                              */
   buffer[buffpos] = '\0';

   /* If haven't used all the buffer, shrink it down to used size       */
   if(++buffpos < bufflen)
   {
      if((buffer = realloc(buffer, bufflen))==NULL)
      {
         errno = ENOMEM;
         return(NULL);
      }
   }
   
   return(buffer);
}
