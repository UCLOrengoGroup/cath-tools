/*************************************************************************

   Program:    
   File:       ParseRes.c
   
   Version:    V1.10R
   Date:       12.10.12
   Function:   Parse a residue specification
   
   Copyright:  (c) SciTech Software 1993-2012
   Author:     Dr. Andrew C. R. Martin
   EMail:      andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
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
   V1.8  29.09.05 Moved ParseResSpec() into DoParseResSpec() with extra
                  param and added wrappers for ParseResSpec() and 
                  ParseResSpecNoUpper()  (Changes by Tony Lewis) By: TL
   V1.9  05.01.12 Default behaviour of ParseResSpec() is now not to
                  upcase the chain label - there are now too many PDB
                  entries with lower case chain names for this to be
                  sensible.   By: ACRM
   V1.10 12.10.12 insert is now a properly terminated string when there is
                  no insert

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
/*>BOOL ParseResSpec(char *spec, char *chain, int *resnum, char *insert)
   ---------------------------------------------------------------------
   Input:   char  *spec    Residue specification
   Output:  char  *chain   Chain label
            int   *resnum  Residue number
            char  *insert  Insert label
   Returns: BOOL           Success?

   Splits up a residue specification of the form 
         [c][.]num[i]
   into chain, resnum and insert. Chain and insert are optional and will
   be set to spaces if not specified. Converts the resiude specification
   to upper case before processing.
   
   Moved the code that was here to a new function, DoParseResSpec()
   and made this function just call that new function.  See
   DoParseResSpec()'s comments for notes on previous changes.  This
   move is to allow the underlying function to have an extra parameter
   to specify whether or not the residue specification should be upper
   cased (without affecting code that calls this function).

   29.09.05 Original   By: TL
   05.01.12 Now behaves the same as ParseResSpecNoUpper(). There are now
            too many PDB files with lower case chain names (e.g. 1gav,
            3n9r, etc.) for the old default behaviour or up-casing 
            everything.   By: ACRM
*/
BOOL ParseResSpec(char *spec, char *chain, int *resnum, char *insert)
{
   return DoParseResSpec(spec, chain, resnum, insert, FALSE);
}

/************************************************************************/
/*>BOOL ParseResSpecNoUpper(char *spec, char *chain, int *resnum, 
                            char *insert)
   --------------------------------------------------------------
   Input:   char  *spec    Residue specification
   Output:  char  *chain   Chain label
            int   *resnum  Residue number
            char  *insert  Insert label
   Returns: BOOL           Success?

   Splits up a residue specification of the form 
         [c][.]num[i]
   into chain, resnum and insert. Chain and insert are optional and will
   be set to spaces if not specified. Does not converts the resiude
   specification to upper case before processing.
   
   29.09.05 Original   By: TL
*/
BOOL ParseResSpecNoUpper(char *spec, char *chain, int *resnum, 
                         char *insert)
{
   return DoParseResSpec(spec, chain, resnum, insert, FALSE);
}

/************************************************************************/
/*>BOOL DoParseResSpec(char *spec, char *chain, int *resnum, char *insert, 
                       BOOL uppercaseresspec)
   -----------------------------------------------------------------------
   Input:   char  *spec    Residue specification
            BOOL           uppercaseresspec
   Output:  char  *chain   Chain label
            int   *resnum  Residue number
            char  *insert  Insert label
   Returns: BOOL           Success?

   Splits up a residue specification of the form 
         [c][.]num[i]
   into chain, resnum and insert. Chain and insert are optional and will
   be set to spaces if not specified. If uppercaseresspec eqauls TRUE,
   the spec is upper cased before processing
   
   21.07.93 Original    By: ACRM
   17.07.95 Added BOOL return
   18.03.98 Added option to include a . to separate chain and residue
            number so numeric chain names can be used
   29.09.05 Moved this code to from ParseResSpec() to DoParseResSpec()
            and made that function just call this new function.
            This move is to allow this underlying function to have an
            extra parameter to specify whether or not the residue
            specification should be upper cased (without affecting code
            that calls the old function). By: TL
   12.10.12 insert is now a properly terminated string when there is
            no insert
*/
BOOL DoParseResSpec(char *spec, char *chain, int *resnum, char *insert, 
                    BOOL uppercaseresspec)
{
   char  *ptr,
         *ptr2;
   BOOL  DoRestore = FALSE,
         retval    = TRUE;

   /* 11.10.99 Default resnum of 0                                      */
   *resnum = 0;

   /* Upper case the residue specification if it has been requested     */
   if (uppercaseresspec == TRUE)
   {
      UPPER(spec);   
   }
   KILLLEADSPACES(ptr, spec);
     
   /* Extract chain from spec                                           */
   if(*ptr == '.')
   {
      *chain = ' ';
      ptr++;
   }
   else if((*(ptr+1) == '.') || (!isdigit(*ptr) && (*ptr != '-')))
   {
      /* Chain was specified                                            */
      *chain = *ptr;
      ptr++;
      if(*ptr == '.')
      {
         ptr++;
      }
   }
   else
   {
      /* Spec started with a digit, so no chain specified               */
      *chain = ' ';
   }
   
   /* Extract insert from spec                                          */
   insert[0] = ' ';
   insert[1] = '\0';  /* Added 12.10.12                                 */
   
   for(ptr2 = ptr; *ptr2; ptr2++)
   {
      /* 11.10.99 Now also checks that it isn't a - as the first 
         character 
      */
      if(!isdigit(*ptr2) && ((ptr2!=ptr)||(*ptr2 != '-')))
      {
         *insert = *ptr2;
         *ptr2   = '\0';
         DoRestore = TRUE;
         break;
      }
   }
   
   /* Extract residue number from spec                                  */
   if(sscanf(ptr,"%d",resnum) == 0)
      retval = FALSE;

   if(DoRestore)
   {
      /* V1.1: Restore the original string                              */
      *ptr2 = *insert;
   }

   return(retval);
}


