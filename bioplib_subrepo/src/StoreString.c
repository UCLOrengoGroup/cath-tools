/*************************************************************************

   Program:    
   File:       StoreString.c
   
   Version:    V1.21
   Date:       18.06.02
   Function:   
   
   Copyright:  (c) SciTech Software 1991-2002
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

*************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include <string.h>
#include "general.h"
#include "macros.h"

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
/*>STRINGLIST *StoreString(STRINGLIST *StringList, char *string)
   -------------------------------------------------------------
   Input:     STRINGLIST  *StringList   The current linked list or NULL 
                                        if nothing yet allocated
              char        *string       The string to store
   Returns:   STRINGLIST  *             Start of linked list. Used on
                                        first call (when input StringList
                                        is NULL) to return the pointer to
                                        the start of the linked list.
                                        NULL if unable to allocate.

   Stores strings (of any length) in a linked list of type STRINGLIST.
   Return a pointer to the start of the linked list which is used on
   the first call to access the newly allocated memory.

   If allocation fails, memory allocated so far is freed and the routine
   returns NULL.

   06.11.95 Original    By: ACRM
*/
STRINGLIST *StoreString(STRINGLIST *StringList, char *string)
{
   STRINGLIST *p, 
              *start;
   
   if((StringList!=NULL) && ((string == NULL) || (string[0] == '\0')))
      return(StringList);
   
   /* Set p to the String List and move to the end of the linked list   */
   start = StringList;

   /* If nothing in the list, initialise it                             */
   if(start == NULL)
   {
      INIT(start,STRINGLIST);
      p=start;
   }
   else  /* Move to end of current list and add another item            */
   {
      p=start;
      LAST(p);

      /* Only allocate another slot if a string is inserted in this one */
      if(p->string != NULL)
         ALLOCNEXT(p,STRINGLIST);
   }
   
   /* Check allocation                                                  */
   if(p==NULL)
   {
      /* If failed, free the list so far and return NULL                */
      FREELIST(start, STRINGLIST);
      return(NULL);
   }
   p->string = NULL;
   
   /* Everything OK, allocate memory for the string                     */
   if((string != NULL) && (string[0] != '\0'))
   {
      if((p->string = (char *)malloc((1+strlen(string))*sizeof(char)))
         ==NULL)
      {
         /* No memory, free linked list and return                      */
         FREELIST(start, STRINGLIST);
         return(NULL);
      }
      
      /* Still OK, copy in the string and return                        */
      strcpy(p->string,string);
   }
   
   return(start);
}


