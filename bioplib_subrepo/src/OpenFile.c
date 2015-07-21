/*************************************************************************

   Program:    
   File:       OpenFile.c
   
   Version:    V1.22
   Date:       28.07.05
   Function:   
   
   Copyright:  (c) SciTech Software 1991-2005
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
   V1.1  08.02.91 Added KillLine()
   V1.2  10.02.91 Added setextn() and index()
   V1.3  20.03.91 Added Word()
   V1.4  28.05.92 ANSIed
   V1.5  22.06.92 Added tab check to Word(). Improved setextn().
                  Added WordN(). Documented other routines.
   V1.6  27.07.93 Corrected fsscanf() for double precision
   V1.7  07.10.93 Checks made on case before toupper()/tolower()
                  for SysV compatibility. Also index() becomes
                  chindex()
   V1.8  18.03.94 getc() -> fgetc()
   V1.9  11.05.94 Added GetFilestem(), upstrcmp(), upstrncmp() &
                  GetWord()
   V1.10 24.08.94 Added OpenStdFiles()
   V1.11 08.03.95 Corrected OpenFile() for non-UNIX
   V1.12 09.03.95 Added check on non-NULL filename in OpenFile()
   V1.13 17.07.95 Added countchar()
   V1.14 18.10.95 Moved YorN() to WindIO.c
   V1.15 06.11.95 Added StoreString(), InStringList() and FreeStringList()
   V1.16 22.11.95 Moved ftostr() to generam.c
   V1.17 15.12.95 Added QueryStrStr()
   V1.18 18.12.95 OpenStdFiles() treats filename of - as stdin/stdout
   V1.19 05.02.96 OpenStdFiles() allows NULL pointers instead if filenames
   V1.20 18.09.96 Added padchar()
   V1.21 18.06.02 Added string.h
   V1.22 28.07.05 Added conditionals for Mac OS/X

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SysDefs.h"
#include "port.h"

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
/*>FILE *OpenFile(char *filename, char *envvar, char *mode, BOOL *noenv)
   ---------------------------------------------------------------------
   Input:     char    *filename     Filename to be opened
              char    *envvar       Unix/MS-DOS environment variable
                                    Other OS assign name (with :)
              char    *mode         Mode in which to open file (r, w, etc)
   Output:    BOOL    *noenv        Set to TRUE under Unix/MS-DOS if 
                                    the reason for failure was that the
                                    environment variable was not set.
   Returns:   FILE    *             File pointer or NULL on failure

   Attempts to open a filename as specified. Returns a file
   pointer. If this fails:

   Under UNIX/MS-DOS:
   gets the contents of the envvar environment variable and prepends
   that to the filename and tries again. If envvar was not set, noenv
   is set to TRUE and the routine returns a NULL pointer.

   Under other OSs:
   prepends the envvar string onto the filename and tries to open the
   file again.

   Returns the pointer returned by the open() command after all this.

   22.09.94 Original    By: ACRM
   11.09.94 Puts a : in for the assign type.
   24.11.94 Added __unix define. Checks for trailing / in environment
            variable
   08.03.95 Corrected basename to filename in non-unix version
   09.03.95 Checks that filename is not a NULL or blank string
   28.07.05 Added conditionals for Mac OS/X: __MACH__ and __APPLE__
*/
FILE *OpenFile(char *filename, char *envvar, char *mode, BOOL *noenv)
{
   char *datadir,
        buffer[160];
   FILE *fp;
   
   if(filename == NULL || filename[0] == '\0')
      return(NULL);

   if(noenv != NULL) *noenv = FALSE;

   /* Try to open the filename as specified                             */
   if((fp=fopen(filename,mode)) == NULL)
   {
      /* Failed, so build alternative directory/filename                */
#if (unix || __unix__ || MS_WINDOWS || __unix || __MACH__ || __APPLE__)
      if((datadir = getenv(envvar)) != NULL)
      {
         if(datadir[strlen(datadir)-1] == '/')
            sprintf(buffer,"%s%s",datadir,filename);
         else
            sprintf(buffer,"%s/%s",datadir,filename);
         fp = fopen(buffer,mode);
      }
      else
      {
         if(noenv != NULL) *noenv = TRUE;
         return(NULL);
      }
#else
      sprintf(buffer,"%s:%s",envvar,filename);
      fp = fopen(buffer,mode);
#endif
   }
   
   return(fp);
}


