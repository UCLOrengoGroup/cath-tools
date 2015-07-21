/*************************************************************************

   Program:    
   File:       NComb.c
   
   Version:    V1.0
   Date:       10.09.96
   Function:   
   
   Copyright:  (c) SciTech Software 1996
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
#include "MathUtil.h"

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
/*>ULONG NComb(int n, int r)
   -------------------------
   Calculates number of combinations of n items in r groups
   Returns 0 if a numeric overflow occurs.
   
   09.09.96 Original   By: ACRM
*/
ULONG NComb(int n, int r)
{
   ULONG f;
   f = factorial(r);
   
   return((f>0L)?(NPerm(n,r)/factorial(r)):0L);
}
