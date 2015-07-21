/*************************************************************************

   Program:    
   File:       IndxReal.c
   
   Version:    V1.3
   Date:       08.07.96
   Function:   Index heapsort a REAL array
   
   Copyright:  (c) SciTech Software 1991-1996
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      andrew@stagleys.demon.co.uk

   Originally written while at:
               Laboratory of Mathematical Biology,
               National Institue for Medical Research,
               The Ridgeway,
               Mill Hill,
               London,
               NW7 1AA
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   This routine uses a heapsort to index a floating point array
   such that arrin[indx[j]] is in ascending order with j.
   It is modified from the FORTRAN version in 'Numerical Recipes'
   Page 233. This version correctly sorts from array element 0
   as opposed to 1 in the FORTRAN version.

**************************************************************************

   Usage:
   ======
   IndexReal(n,arrin,indx)
   Input:  int   n       Number of elements in array
           REAL  *arrin  Array to be indexed
   Output: int   *indx   Index array

**************************************************************************

   Revision History:
   =================
   V1.0  23.06.90 Original
   V1.1  01.06.92 ANSIed and autodoc'd
   V1.2  19.07.93 Corrected bug (said j=l+1 rather than j=l+l). Oops!
   V1.3  08.07.96 Changed from double to REAL and tidied.

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
/*>void IndexReal(REAL *arrin, int *indx, int n)
   ---------------------------------------------
   Input:  REAL   *arrin   Array to be indexed
           int    n        Number of elements in array
   Output: int    *indx    Index array
   
   Index an array by Heapsort.
   
   03.06.90 Original
   01.06.92 ANSIed and autodoc'd
   19.07.93 Corrected bug (said j=l+1 rather than j=l+l). Oops!
   08.07.96 Changed from double to REAL and tidied. Changed param order
*/
void IndexReal(REAL *arrin, int *indx, int n)
{
   int  i, j, l,
        ir,
        indxt;
   REAL q;

   for(j=0; j<n; j++) 
      indx[j]=j;

   l  = n/2+1;
   ir = n;

   for(;;)
   {
      if(l>1)
      {
         indxt = indx[--l - 1];
         q     = arrin[indxt];
      }
      else
      {
         indxt      = indx[ir-1];
         q          = arrin[indxt];
         indx[ir-1] = indx[0];
         if(--ir == 1)
         {
            indx[0] = indxt;
            return;
         }
      }
      i = l;
      j = l+l;

      while(j <= ir)
      {
         if(j < ir)
         {
            if(arrin[indx[j-1]] < arrin[indx[j]]) 
               j++;
         }

         if(q < arrin[indx[j-1]])
         {
            indx[i-1] = indx[j-1];
            i         = j;
            j        += j;
         }
         else
         {
            j = ir+1;
         }
      }
      indx[i-1] = indxt;
   }
}

