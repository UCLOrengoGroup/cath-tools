/*************************************************************************

   Program:    
   File:       strcatalloc.c
   
   Version:    V1.1
   Date:       11.07.00
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
   V1.0  22.05.99 Original   By: ACRM
   V1.1  11.07.00 Check that realloc succeeded

*************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include <string.h>

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
/*>char *strcatalloc(char *instr, char *catstr)
   --------------------------------------------
   Input:   char    *instr    String to append to
            char    *catstr   String to append
   Returns: char    *         realloc'd version of instr with catstr
                              appended

   Like strcat() but uses a realloc() on instr to make space available.

   22.05.99 Original   By: ACRM
   16.06.99 Initialise outstr to NULL
   25.08.99 Fixed bug where testing for NULL outstr instead of catstr
   11.07.00 Check that realloc succeeded
*/
char *strcatalloc(char *instr, char *catstr)
{
   int  totLen;
   char *outstr = NULL;
   
   totLen = ((instr==NULL)  ? 0 : strlen(instr)) + 
      ((catstr==NULL) ? 0 : strlen(catstr));
   if((outstr = realloc(instr, totLen+1))!=NULL)
   {
      if(instr==NULL)
         outstr[0] = '\0';
      strcat(outstr, catstr);
   }
   
   return(outstr);
}
