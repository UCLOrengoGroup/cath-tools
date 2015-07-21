/*************************************************************************

   Program:    
   File:       justify.c
   
   Version:    V1.1
   Date:       18.06.02
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
   V1.0  30.05.02 Original
   V1.1  18.06.02 Added string.h

*************************************************************************/
/* Includes
*/
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
/*>void RightJustify(char *string)
   -------------------------------
   I/O:     char    *string           A string 

   Right justifies a string in place

   30.05.02 Original   By: ACRM
*/
void RightJustify(char *string)
{
   int len, dec;

   len = strlen(string);
   len--;
   
   if(len)
   {
      while(string[len] == ' ')
      {
         for(dec = len; dec; dec--)
         {
            string[dec] = string[dec-1];
         }
         string[0] = ' ';
      }
   }
}

