/*************************************************************************

   Program:    
   File:       pearson1.c
   
   Version:    V1.2
   Date:       06.10.98
   Function:   
   
   Copyright:  (c) Dr. Andrew C. R. Martin, UCL, 1996-8
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
/*>REAL pearson1(REAL *x, REAL *y, int NItem)
   ------------------------------------------
   Input:   REAL  *x     Array of x items
            REAL  *x     Array of y items
            int   NItem  Number of items
   Returns: REAL         Pearson correlation coefficient

   This version makes a single pass through the data

   15.07.94 Original    By: ACRM
   18.01.96 Alternative version which does a single pass through the
            data. Method from page 191 of Statistical Methods in 
            Biology 2 ed. Norman TJ Bailey. Publ. Unibooks 1981  By: RM
*/
REAL pearson1(REAL *x, REAL *y, int NItem)
{
   REAL SumX,
        SumY,
        SumXSq,
        SumYSq,
        SumXY,
        Sx, Sy, c, n,
        r;
   int  i;

   SumX = SumY = SumXSq = SumYSq = SumXY = (REAL)0.0;
   for(i=0; i<NItem; i++)
   {
      SumX   += x[i];
      SumY   += y[i];
      SumXSq += x[i]*x[i];
      SumYSq += y[i]*y[i];
      SumXY  += x[i]*y[i];
   }
   /* Calculate correlation coefficient                                 */
   
   n = (REAL)NItem;
   
   Sx = SumXSq - SumX*SumX/n;
   Sy = SumYSq - SumY*SumY/n;
   c  = SumXY  - SumX*SumY/n;
   
   r = c / (REAL)sqrt((double)(Sx * Sy));
   
   return(r);
}


