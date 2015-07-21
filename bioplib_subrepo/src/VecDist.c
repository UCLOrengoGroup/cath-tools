/*************************************************************************

   Program:    
   File:       VecDist.c
   
   Version:    V1.2
   Date:       06.10.98
   Function:   
   
   Copyright:  (c) SciTech Software 1996-8
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
   V1.0  29.01.96 Original   By: ACRM
   V1.1  18.06.96 Added vector routines
   V1.2  06.10.98 Added VecAdd3()

*************************************************************************/
/* Includes
*/
#include <math.h>
#include "MathType.h"

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
/*>REAL VecDist(REAL *a, REAL *b, int len)
   ---------------------------------------
   Input:   REAL    *a     An arbitrary length vector (as an array)
            REAL    *b     An arbitrary length vector (as an array)
            int     len    The dimensionality of the vectors (array
                           length)
   Returns: REAL           The distance between the points described by
                           the two vectors

   Finds the distance between two vectors of arbitrary length

   28.07.95 Original    By: ACRM
*/
REAL VecDist(REAL *a, REAL *b, int len)
{
   REAL sumsq = 0.0;
   int  i;
   
   for(i=0; i<len; i++)
      sumsq += (a[i] - b[i]) * (a[i] - b[i]);

   return((REAL)sqrt((double)sumsq));
}


