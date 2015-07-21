/*************************************************************************

   Program:    
   File:       KillLeadSpaces.c
   
   Version:    V1.20
   Date:       18.09.96
   Function:   
   
   Copyright:  (c) SciTech Software 1991-6
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

*************************************************************************/
/* Includes
*/

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
/*>char *KillLeadSpaces(char *string)
   ----------------------------------
   Input:      char  *string        A character string
   Returns:   *char                 A pointer to the string with the 
                                    leading spaces removed

   This routine strips leading spaces and tabs from a string returning
   a pointer to the first non-whitespace character.

   N.B. THE MACRO KILLLEADSPACES() MAY NOW BE USED INSTEAD
   
   06.02.91 Original
   28.05.92 ANSIed
   06.07.93 Added tab skipping
*/
char *KillLeadSpaces(char *string)
{
   while (*string == ' ' || *string == '\t') string++;
   return(string);
}


