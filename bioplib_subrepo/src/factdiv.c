/*************************************************************************

   Program:    
   File:       factdiv.c
   
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
/*>ULONG factdiv(int n1, int n2)
   -----------------------------
   Calculates the factorial of one number divided by the factorial of
   another (smaller) number.
   Returns 0 on numeric overflow or if n2 > n1

   09.09.96 Original   By: ACRM
*/
ULONG factdiv(int n1, int n2)
{
   int   i;
   ULONG ret  = 1L,
         prev = 0L;

   if(n2 > n1)
      return(0);
   
   for(i=n2+1; i<=n1; i++)
   {
      ret *= (ULONG)i;
      if(ret < prev)
         return(0);
      prev = ret;
   }
   
   return(ret);
}


