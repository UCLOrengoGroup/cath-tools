/*************************************************************************

   Program:    
   File:       getfield.c
   
   Version:    V1.1
   Date:       18.06.02
   Function:   
   
   Copyright:  (c) Dr. Andrew C. R. Martin, University of Reading, 2002
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
/*>void getfield(char *buffer, int start, int width, char *str)
   ------------------------------------------------------------
   Input:     char   *buffer      Buffer from which to read a field
              int    start        Starting column (count from 0)
              int    width        Width of field to read
   Output:    char   *str         Field read from buffer

   Reads a column out of a buffer. If the specfied column extends beyond
   the size of the buffer, then it will be padded with spaces.

   Note that the output string must be of at lease width+1 characters
   to store the field read from the buffer plus the terminating
   character.

   30.05.02 Original   By: ACRM
*/
void getfield(char *buffer, int start, int width, char *str)
{
   int i, 
       j,
       len;
   
   len = strlen(buffer);
   
   for(i=0, j=0; i<width; i++)
   {
      if(start+i >= len)
      {
         str[j++] = ' ';
      }
      else
      {
         str[j++] = buffer[start+i];
      }
   }
   str[j] = '\0';
}

   
