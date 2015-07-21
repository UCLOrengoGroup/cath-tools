/*************************************************************************

   Program:    
   File:       CheckExtn.c
   
   Version:    V1.20R
   Date:       18.09.96
   Function:   General purpose routines
   
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
   These are generally useful C routines, mostly string handling.

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

*************************************************************************/
/* Includes
*/

#include <string.h>
#include <stdlib.h>
#include "general.h"
#include "SysDefs.h"


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
/*>BOOL CheckExtn(char *string, char *ext)
   ---------------------------------------
   Input:   char *string    String to be checked for given extension
            char *extn      Extension to check for
   Returns: BOOL            Found?

   Check the extension of a filename. For use on machines like VAXes,
   MS-DOS and Amigas, everything is converted to upper case first.
   18.06.93 Original    By: ACRM
*/
BOOL CheckExtn(char  *string,
               char  *ext)
{
   int   extl     = strlen(ext),
         strl     = strlen(string);
   char  *buff1,
         *buff2;
   BOOL  RetVal   = TRUE;
         
   buff1 = (char *)malloc(strl * sizeof(char));
   buff2 = (char *)malloc(extl * sizeof(char));
   
   if(buff1==NULL || buff2==NULL)
   {
      if(buff1) free(buff1);
      if(buff2) free(buff2);
      return(FALSE);
   }
         
   StringToUpper(string,buff1);
   StringToUpper(ext,   buff2);
         
   if(strncmp(buff2,buff1+(strl-extl),extl))
      RetVal = FALSE;
      
   free(buff1);
   free(buff2);
   return(RetVal);
}


