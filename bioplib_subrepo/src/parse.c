/*************************************************************************

   Program:    
   File:       parse.c
   
   Version:    V1.10R
   Date:       28.02.11
   Function:   A keyword command parser
   
   Copyright:  (c) SciTech Software 1990-2011
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
   parse() is a command line parser which will accept upper or
   lower case commands and abbreviations. Comment lines may be
   indicated using a !. The keyword structure array and returned
   string array are defined thus:
           KeyWd keywords[NCOMM];
           char  *strparam[MAXSTRPARAM];
   The returned REAL parameters are defined thus:
           REAL  floatparam[MAXFLOATPARAM];
   Space for the returned strings must be allocated thus:
           strparam[n] = (char *)malloc(MAXSTRLEN * sizeof(char));
   and repeated for each parameter.

   The keyword list with type and numbers of returned parameters
   is constructed using the MAKEKEY macros:
           MAKEKEY(keywords[0],"RANGE",NUMBER,2);
           MAKEKEY(keywords[1],"STRING",STRING,1);
   Here, the keywords must be defined in upper case.


   mparse() is used in the same way, but allows a variable number of
   parameters for each keyword. Keywords are of type MKeyWd and are
   defined using the macro MAKEMKEY:
           MAKEMKEY(keywords[0],"RANGE",NUMBER,2,2);
           MAKEMKEY(keywords[1],"STRING",STRING,1,3);

**************************************************************************

   Usage:
   ======
   parse(comline,nkeys,keywords,floatparam,strparam)
   Input:      char  *comline       A command line string to parse
               int   nkeys          Number of keywords
               KeyWd *keywords      Array of keyword structures
   Output:     REAL  *floatparam    Array of returned strings
               char  **strparam     Array of pointers to returned strings
   Returns:    int                  Index of found command or error flag

   mparse(comline,nkeys,keywords,floatparam,strparam,nparam)
   Input:      char   *comline       A command line string to parse
               int    nkeys          Number of keywords
               MKeyWd *keywords      Array of keyword structures
   Output:     REAL   *floatparam    Array of returned strings
               char   **strparam     Array of pointers to returned strings
               int    *nparam        Number of parameters found
   Returns:    int                   Index of found command or error flag


**************************************************************************

   Revision History:
   =================
   V1.0  11.07.90 Original
   V1.1  29.10.90 match() now frees the memory it allocates and calls
                  terminate()
                  parse() now calls terminate() on the keyword string
   V1.2  25.09.91 Messages will only appear from parse() if NOISY is
                  #defined.
                  Added FPU support.
   V1.3  28.05.92 ANSIed and autodoc'd
   V1.4  08.12.92 Includes stdlib.h
   V1.5  22.04.93 Various tidying to exact ANSI standard and of function
                  headers. Corrected some calls to free()
   V1.6  16.06.93 Tidied for book
   V1.7  01.03.94 Added mparse()
   V1.8  11.03.94 Added internal support for lines starting with a $.
                  The line is passed as a system() call and parse()
                  acts as if the line had been a comment.
   V1.9  08.10.99 Initialised some variables
   V1.10 28.02.11 Added # as a comment introducer

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

#include "MathType.h"
#include "SysDefs.h"
#include "macros.h"
#include "parse.h"
#include "general.h"

/************************************************************************/
/* General defines for these routines
*/
#define LF  10
#define CR  13
#define DIC 34       /* Double inverted commas                          */

/************************************************************************/
/*>int parse(char *comline, int nkeys, KeyWd *keywords,
             REAL *floatparam, char **strparam)
   ----------------------------------------------------
   Input:      char  *comline       A command line string to parse
               int   nkeys          Number of keywords
               KeyWd *keywords      Array of keyword structures
   Output:     REAL  *floatparam    Array of returned strings
               char  **strparam     Array of pointers to returned strings
   Returns:    int                  Index of found command or error flag

   Keyword-based command parser. Fixed number of parameters.

   11.07.90 Original    By: ACRM
   22.04.93 Tidied comments, etc. Corrected NULL to 0.
   11.03.94 Added $ line handling
   08.10.99 Initialise nlett
   28.02.11 Added # for comments
*/
int parse(char  *comline,
          int   nkeys,
          KeyWd *keywords,
          REAL  *floatparam,
          char  **strparam)
{
   char *command;
   int  i,n,found,nletters,nlett = 0;

   command = KillLeadSpaces(comline);
   TERMINATE(command);

   if(command[0] == '$')
   {
      system(command+1);
      return(PARSE_COMMENT);
   }
   
   found = 0;
   if((command[0]=='!')  ||
      (command[0]=='#')  ||
      (command[0]==LF)   ||
      (command[0]==CR)   ||
      (command[0]=='\0'))
      return(PARSE_COMMENT);

   for(i=0;i<nkeys;i++)
   {
      /* match() returns 1 if first string finishes first or exact match
                         2 if second string finishes first
                         0 if a mismatch
         We only want to act in the first case
      */
      if((n=match(command,(keywords[i]).name,&nletters))==1)
      {
         if(found)      /* If found already                             */
         {
            return(PARSE_ERRC);
         }
         found = i+1;   /* +1, so keyword 0 will flag TRUE              */
         nlett = nletters;
      }
   }
   if(!found)
   {
      return(PARSE_ERRC);
   }
   command+=nlett;
   found--;             /* Reset to point to the correct keyword        */

   /* Get data requirements for this keyword                            */
   if((keywords[found]).string)
   {
      for(i=0; i<(keywords[found]).nparam; i++)
      {
         command = KillLeadSpaces(command);
         if((nletters = GetString(command,strparam[i]))==0)
         {
            return(PARSE_ERRP);
         }
         command += nletters;
      }  /* End of for(i)                                               */
   }
   else
   {
      /* A numeric or no parameter                                      */
      for(i=0; i<(keywords[found]).nparam; i++)
      {
         command = KillLeadSpaces(command);
         if(!GetParam(command,&(floatparam[i]),&nletters))
         {
            return(PARSE_ERRP);
         }
         command += nletters;
      }  /* End of for(i)                                               */
   }  /* End of else                                                    */
   return(found);
}

/************************************************************************/
/*>int match(char *comstring, char *string2, int *nletters)
   --------------------------------------------------------
   Input:      char  *comstring     A character string
               char  *string2       A second string
   Output:     int   *nletters      Number of letters matched
   Returns:    int                  0 String mismatch
                                    1 First string finished first
                                    2 Second string finished first

   This routine matches two strings, but stops the comparison as soon
   as a space or NULL is found in either string. The returned value
   indicates which string finished first or 0 if the letters before the
   space or NULL have a mismatch. The routine calls StringToUpper()
   on `comstring' before the comparison.

   11.07.90 Original    By: ACRM
   22.04.93 Tidied comments, etc. Added check on malloc and corrected
            calls to free()
*/
int match(char *comstring,
          char *string2,
          int  *nletters)
{
   int  i;
   char *string1;

   TERMINATE(comstring);
   TERMINATE(string2);
   string1 = (char *)malloc((strlen(comstring) + 2) * sizeof(char));
   if(string1 == NULL) return(0);

   StringToUpper(comstring,string1);

   for(i=0;;i++)
   {
      if((!string1[i])||(string1[i]==' '))
      {
         *nletters = i;
         free(string1);
         return(1);
      }
      if((!string2[i])||(string2[i]==' '))
      {
         *nletters = i;
         free(string1);
         return(2);
      }
      if(string1[i] != string2[i])
      {
         *nletters = i;
         free(string1);
         return(0);
      }
   }
}

/************************************************************************/
/*>int GetString(char *command, char *strparam)
   --------------------------------------------
   Input:      char  *command       A character string
   Output:     char  *strparam      Returned character string
   Returns:    int                  Number of characters pulled out
                                    of the command string

   This routine returns the first space-delimited group of characters
   from character string `command'

   11.07.90 Original    By: ACRM
   22.04.93 Tidied comments, etc. Changed toggle method
*/
int GetString(char *command,
              char *strparam)
{
   int i,j,inv_commas;

   inv_commas=0;
   j=0;
   for(i=0;;i++)
   {
      if(command[i]==DIC)
      {
         /* Toggle the inv_commas flag                                  */
         inv_commas = !inv_commas;

         /* Don't copy anything                                         */
         continue;
      }

      /* Break out if we're at the end of a line                        */
      if((command[i]==LF)
       ||(command[i]==CR)
       ||(command[i]=='\0')) break;

      /* Also break out if we've a space and we're not between
         inverted commas
      */
      if((command[i]==' ') && (!inv_commas)) break;

      /* Other wise copy the character                                  */
      strparam[j++] = command[i];
   }
   strparam[j]='\0';
   return(i);
}

/************************************************************************/
/*>int GetParam(char *command, REAL *value, int *nletters)
   -------------------------------------------------------
   Input:      char  *command       A character string
   Output:     REAL  *value         Returned float value
               int   *nletters      Number of charcters pulled out
                                    of the command string
   Returns:    int                  0 If error
                                    1 If OK

   This routine extracts the first space-delimited number from the
   `command' character string.

   11.07.90 Original    By: ACRM
   22.04.93 Tidied comments, etc. Corrected NULL to 0
*/
int GetParam(char  *command,
             REAL  *value,
             int   *nletters)
{
   char buffer[50];
   int  retval;

   if((*nletters = GetString(command,buffer))==0)
      return(0);

   retval = sscanf(buffer,"%lf",value);
   return(retval);
}

/************************************************************************/
/*>int mparse(char *comline, int nkeys, MKeyWd *keywords,
              REAL *floatparam, char **strparam, int *nparam)
   ----------------------------------------------------------
   Input:      char   *comline       A command line string to parse
               int    nkeys          Number of keywords
               MKeyWd *keywords      Array of keyword structures
   Output:     REAL   *floatparam    Array of returned strings
               char   **strparam     Array of pointers to returned strings
               int    *nparam        Number of parameters found
   Returns:    int                   Index of found command or error flag

   As parse(), but allows variable number of parameters to each keyword.

   23.02.94 Original based on parse()   By: ACRM
   11.03.94 Added $ line handling
   08.10.99 Initialise nlett to 0
*/
int mparse(char   *comline,
           int    nkeys,
           MKeyWd *keywords,
           REAL   *floatparam,
           char   **strparam,
           int    *nparam)
{
   char *command;
   int  i,n,found,nletters,nlett=0;

   command = KillLeadSpaces(comline);
   TERMINATE(command);
   
   if(command[0] == '$')
   {
      system(command+1);
      return(PARSE_COMMENT);
   }

   found = 0;
   if((command[0]=='!') ||
      (command[0]==LF)  ||
      (command[0]==CR)  ||
      (command[0]=='\0'))
      return(PARSE_COMMENT);

   for(i=0;i<nkeys;i++)
   {
      /* match() returns 1 if first string finishes first or exact match
                         2 if second string finishes first
                         0 if a mismatch
         We only want to act in the first case
      */
      if((n=match(command,(keywords[i]).name,&nletters))==1)
      {
         if(found)      /* If found already                             */
         {
            return(PARSE_ERRC);
         }
         found = i+1;   /* +1, so keyword 0 will flag TRUE              */
         nlett = nletters;
      }
   }

   if(!found)
   {
      return(PARSE_ERRC);
   }

   command+=nlett;
   found--;             /* Reset to point to the correct keyword        */

   *nparam = 0;         /* Zero the parameter count                     */

   /* Get data requirements for this keyword                            */
   if((keywords[found]).string)
   {
      for(i=0; i<(keywords[found]).maxparam; i++)
      {
         command = KillLeadSpaces(command);
         if((nletters = GetString(command,strparam[i]))==0)
         {
            if(i < (keywords[found]).minparam)
               return(PARSE_ERRP);
            else
               break;
         }
         else
         {
            (*nparam)++;
         }
         command += nletters;
      }  /* End of for(i)                                               */
   }
   else
   {
      /* A numeric or no parameter                                      */
      for(i=0; i<(keywords[found]).maxparam; i++)
      {
         command = KillLeadSpaces(command);
         if(!GetParam(command,&(floatparam[i]),&nletters))
         {
            if(i < (keywords[found]).minparam)
               return(PARSE_ERRP);
            else
               break;
         }
         command += nletters;
         (*nparam)++;
      }  /* End of for(i)                                               */
   }  /* End of else                                                    */
   return(found);
}

