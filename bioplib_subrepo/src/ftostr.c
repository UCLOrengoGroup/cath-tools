/*************************************************************************

   Program:    
   File:       ftostr.c
   
   Version:    V1.3R
   Date:       03.06.05
   Function:   Convert a REAL to a string
   
   Copyright:  (c) SciTech Software 1991-2005
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
   Convert a REAL to a string

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  22.11.95 Original - Routines from general.c
   V1.1  09.10.97 Fixed inconsistent handling of %.0f between C libs
   V1.2  14.11.97 Another fix to same
   V1.3  03.06.05 Tidied up to stop warnings under GCC 3.2.2

*************************************************************************/
/* Includes
*/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MathType.h"
#include "SysDefs.h"
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
/*>char *ftostr(char *str, int maxlen, REAL x, int precision)
   ----------------------------------------------------------
   Output:  char  *str       String version of number
   Input:   int   maxlen     Max length for f-form. Switches to e-form
                             if exceeded
            REAL  x          Number to convert
            int   precision  Number of decimal places
                             If negative, use e-form rather than f-form
   Returns: char  *          Same as str

   Convert a REAL to a string using precision decimal places. If 
   precision is negative, use e-form, otherwise use f-form. This is used
   to generate precisely formatted string versions of numbers for
   applications where the appearance of a numeric value is important.

   01.07.92 Original
   07.07.92 Changed to stop value of -0.0
   24.07.92 Added maxlen parameter. If the f-form will exceed this number
            of characters we switch to e-form autmatically
   28.07.92 If precision -ve, is never < maxlen-5
   27.07.93 Changed to double precision I/O
   09.10.97 If precision is 0 prints as an integer since C libraries
            handle "%.0f" inconsistently
   14.11.97 Oops, had forgotten to fix this in the check for e-form
   03.06.05 Tidied up some loops to stop warnings under GCC 3.2.2
*/
char *ftostr(char  *str,
             int   maxlen,
             REAL  x,
             int   precision)
{
   char   fmt[8],
          *ptr;
   int    i;
   REAL   val;
   
   /* If usiing e-form, check whether string will fit                   */
   if(precision > 0)
   {
      REAL lgx;
      int  len;
      
      lgx = log10(ABS(x));       /* Number of pre-dp places - 1         */
      if(lgx < 0.0) lgx = 0.1;   /* i.e. If x<1.0, 1 p-dp after rounding*/
      len = (int)ceil(lgx);      /* Round up to number of pre-dp places */
      len += precision + 2;      /* dp and '\0'                         */
      
      if(len > maxlen)           /* String too long                     */
      {
         if(precision == 0)   precision  = -1;  /* Set to -1            */
         else                 precision *= -1;  /* Switch to e-form     */
      }
   }
   
   /* 09.10.97 Added check on 0 precision                               */
   if((precision != 0) && (precision < -1*(maxlen-5))) 
      precision = -1*(maxlen-5);

   if(precision > 0)             /* 09.10.97 Changed from >=            */
   {
      sprintf(fmt,"%%.%dlf",precision);
   }
   else
   {
      precision = ABS(precision);
      sprintf(fmt,"%%.%dle",precision);
   }
   if(precision == 0)            /* Added 09.10.97                      */
   {
      sprintf(str,"%d",(int)(x+0.5));
   }
   else
   {
      sprintf(str,fmt,x);
   }
   

   if(!strchr(str,'e'))                   /* Not e-form                 */
   {
      if(precision == 0)
      {
         /* Terminate string at .                                       */
         ptr = str;
         while(*ptr)
         {
            if(*ptr++ == '.')
            {
               *ptr = '\0';
               break;
            }
         }
      }
      else
      {
         if((ptr = strchr(str,'.'))!=NULL)   /* Contains a decimal point*/
         {
            ptr++;
            for(i=0; *ptr; ptr++, i++);
            if(i<precision)               /* Fewer places than precision*/
            {
               /* Pad with zeros                                        */
               for( ; i<precision; i++) *ptr++ = '0';
               *ptr = '\0';
            }
         }
         else                             /* No decimal point           */
         {
            /* Move to end of string                                    */
            for(ptr = str; *ptr; ptr++);
            /* Put in a dp                                              */
            *ptr++ = '.';
            /* Pad with zeros                                           */
            for(i=0; i<precision; i++) *ptr++ = '0';
            *ptr = '\0';
         }
      }
   }
      
   if(str[strlen(str)-1] == '.') str[strlen(str)-1] = '\0';
   
   /* Strip any leading spaces                                          */
   /* 03.06.05 Tidied loop for gcc 3.2.2                                */
   /* 30.09.05 Fixed the warning: operation on `i' may be undefined     */
   while(str[0] == ' ')
   {
      i=0;
      while(i<strlen(str))
      {
         str[i] = str[i+1];
         i++;
      }
   }
   

   /* Remove any leading minus sign if value is 0.0                     */
   sscanf(str,"%lf",&val);
   if(val == 0.0 && str[0] == '-')
   {
      /* 03.06.05 Tidied loop for gcc 3.2.2                             */
      /* 30.09.05 Fixed the warning: operation on `i' may be undefined  */
      i=0;
      while(i<strlen(str))
      {
         str[i] = str[i+1];
         i++;
      }
   }

   return(str);
}

