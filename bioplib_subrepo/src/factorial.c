/*************************************************************************

   Program:    
   File:       factorial.c
   
   Version:    V1.0
   Date:       10.09.96
   Function:   
   
   Copyright:  (c) Dr. Andrew C. R. Martin, UCL, 1996
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
   V1.0  10.09.96 Original

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
/*>ULONG factorial(int n)
   ----------------------
   Calculates the factorial of an integer.
   Returns 0 on numeric overflow.
   
   09.09.96 Original   By: ACRM
*/
ULONG factorial(int n)
{
   int i;
   ULONG ret  = 1L,
         prev = 0L;
   
   for(i=2; i<=n; i++)
   {
      ret *= (ULONG)i;
      if(ret < prev)
         return(0);
      prev = ret;
   }

   return(ret);
}


