/*************************************************************************

   Program:    
   File:       AtomNameMatch.c
   
   Version:    V1.7R
   Date:       11.10.99
   Function:   Tests for matching atom names with wild cards
   
   Copyright:  (c) SciTech Software 1993-9
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
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
   V1.0  01.03.94 Original
   V1.1  07.07.95 Now non-destructive
   V1.2  17.07.95 Now checks that a number was specified as part of the
                  spec. and returns a BOOL
   V1.3  23.10.95 Moved FindResidueSpec() from PDBList.c
   V1.4  08.02.96 Added FindResidue() and changed FindResidueSpec() to
                  use it
   V1.5  23.07.96 Added AtomNameMatch() and LegalAtomSpec()
   V1.6  18.03.98 Added option to include a . to separate chain and 
                  residue number so numeric chain names can be used
   V1.7  11.10.99 Allow a . to be used to start a number (such that the
                  default blank chain name is used). Allows negative 
                  residue numbers

*************************************************************************/
/* Includes
*/
#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "SysDefs.h"
#include "pdb.h"

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
/*>BOOL AtomNameMatch(char *atnam, char *spec, BOOL *ErrorWarn)
   ------------------------------------------------------------
   Input:   char   *atnam      The atom name to test
            char   *spec       The atom specification
   I/O:     BOOL   *ErrorWarn  On input, if TRUE, this routine will
                               indicate errors.
                               On output, indicates whether there
                               was an error.
                               Note that you must be careful to supply
                               an lvalue here, you can't just use TRUE
                               or FALSE since it's modified on return.
                               NULL is allowed if you don't care about
                               errors.

   Tests whether an atom name matches an atom name specification.
   ? or % is used to match a single character
   * is used to match any trailing characters; it may not be used for
   leading characters or in the middle of a specification (e.g. *B*,
   C*2 are both illegal).
   Wildcards may be escaped with a backslash.

   For example: C* matches all carbon atoms,
                O5\* matches an atom called O5*
                ?B* matches all beta atoms

   23.07.96 Original   By: ACRM
*/
BOOL AtomNameMatch(char *atnam, char *spec, BOOL *ErrorWarn)
{
   char *specp,
        *atnamp;
   
   /* Step through the specification and the atom name                  */
   for(specp=spec, atnamp = atnam; *specp; specp++, atnamp++)
   {
      switch(*specp)
      {
      case '\\':
         /* If the specification has a \ then we are escaping the next
            character, so just step on to that character
         */
         specp++;
         break;
      case '?':
         /* A query in the specification matches anything, so just
            continue
         */
         continue;
      case '*':
         /* Matches the rest of the string                              */
         if(ErrorWarn != NULL)
         {
            /* Check that there aren't any illegal characters following */
            if(*(specp+1) && *(specp+1) != ' ')
            {
               if(*ErrorWarn)
               {
                  fprintf(stderr,"Error in atom wildcard: %s\n",spec);
               }
               *ErrorWarn = TRUE;
            }
            else
            {
               *ErrorWarn = FALSE;
            }
         }
         return(TRUE);
      default:
         break;
      }

      /* If there is a mismatch return FALSE                            */
      if(*specp != *atnamp)
      {
         if(ErrorWarn != NULL)
            *ErrorWarn = FALSE;
         return(FALSE);
      }

      /* 07.06.05 If both specifications have ended with a space of 
         end of string then return TRUE. Fixed for if the atnam is
         shorter (after moving the alternate atom indicator into its
         own field)
      */
      if((*specp == ' ') && ((*atnamp == ' ') || (*atnamp == '\0')))
      {
         if(ErrorWarn != NULL)
            *ErrorWarn = FALSE;
         return(TRUE);
      }
   }

   /* There have been no errors and we don't need the error flag again  */
   if(ErrorWarn != NULL)
      *ErrorWarn = FALSE;

   /* The specification has run out, see if there are any atom characters
      left
   */
   if(*atnamp && *atnamp!=' ')
      return(FALSE);

   /* Both have ended OK, so the names match                            */
   return(TRUE);
}


/************************************************************************/
/*>BOOL AtomNameRawMatch(char *atnam, char *spec, BOOL *ErrorWarn)
   ---------------------------------------------------------------
   Input:   char   *atnam      The atom name to check
            char   *spec       The atom specification
   I/O:     BOOL   *ErrorWarn  On input, if TRUE, this routine will
                               indicate errors.
                               On output, indicates whether there
                               was an error.
                               Note that you must be careful to supply
                               an lvalue here, you can't just use TRUE
                               or FALSE since it's modified on return.
                               NULL is allowed if you don't care about
                               errors.

   Tests whether an atom name matches an atom name specification.

   This version should be given the raw atom name rather than the 
   massaged one. i.e. " CA " is C-alpha, "CA  " is Calcium

   Normally it checks against the second character onwards unless the
   spec starts with a < in which case it checks from the beginning of
   the string

   Written as a wrapper to AtomNameMatch()

   15.02.01 Original   By: ACRM
*/
BOOL AtomNameRawMatch(char *atnam, char *spec, BOOL *ErrorWarn)
{
   /* If atom spec starts with a < then just bump the spec pointer, 
      otherwise bump the atom name pointer since we will look from the 
      second character of the atom name
   */
   if(*spec == '<')
   {
      spec++;
   }
   else
   {
      atnam++;
   }

   return(AtomNameMatch(atnam, spec, ErrorWarn));
}

#ifdef TEST_MAIN
int main(int argc, char **argv)
{
   char spec[8], atnam[8];
   
   strcpy(atnam, " CA*");
   printf("Atom name '%s':\n", atnam);

   strcpy(spec,"CA");
   printf("'%s' matches? %s\n", spec, (AtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));
   
   strcpy(spec,"<CA");
   printf("'%s' matches? %s\n", spec, (AtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));
   
   strcpy(spec,"C*");
   printf("'%s' matches? %s\n", spec, (AtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));
   
   strcpy(spec,"CA*");
   printf("'%s' matches? %s\n", spec, (AtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));
   
   strcpy(spec,"CA?");
   printf("'%s' matches? %s\n", spec, (AtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));
   
   strcpy(spec,"C\\*");
   printf("'%s' matches? %s\n", spec, (AtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));
   
   strcpy(spec,"C?");
   printf("'%s' matches? %s\n", spec, (AtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));

   return(0);
}
#endif
