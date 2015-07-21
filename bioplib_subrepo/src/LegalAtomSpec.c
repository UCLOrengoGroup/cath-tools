/*************************************************************************

   Program:    
   File:       LegalAtomSpec.c
   
   Version:    V1.7
   Date:       11.10.99
   Function:   
   
   Copyright:  (c) Dr. Andrew C. R. Martin, University of Reading, 2002
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
/*>BOOL LegalAtomSpec(char *spec)
   ------------------------------
   Partner routine for AtomNameMatch(). Checks whether a wildcard
   specfication is legal (i.e. will not return an error when used
   with AtomNameMatch()).

   The only thing which is not legal is characters following a *

   23.07.96 Original   By: ACRM
*/
BOOL LegalAtomSpec(char *spec)
{
   char *chp;
   
   for(chp=spec; *chp; chp++)
   {
      if(*chp == '\\')
      {
         chp++;
      }
      else if(*chp == '*')
      {
         chp++;
         if(*chp && *chp != ' ')
            return(FALSE);
      }
   }
   return(TRUE);
}

