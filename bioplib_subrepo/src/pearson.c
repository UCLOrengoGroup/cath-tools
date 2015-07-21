/*************************************************************************

   Program:    
   File:       pearson.c
   
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
/*>REAL pearson(REAL *x, REAL *y, int NItem)
   -----------------------------------------
   Input:   REAL  *x     Array of x items
            REAL  *x     Array of y items
            int   NItem  Number of items
   Returns: REAL         Pearson correlation coefficient

   This version makes 2 passes through the data

   15.07.94 Original    By: ACRM
*/
REAL pearson(REAL *x, REAL *y, int NItem)
{
   REAL MeanX,
        MeanY,
        SumXSq,
        SumYSq,
        SumXY,
        xt, yt,
        r;
   int  i;

   /* Calculate means of x and y                                        */
   MeanX = MeanY = (REAL)0.0;
   for(i=0; i<NItem; i++)
   {
      MeanX += x[i];
      MeanY += y[i];
   }
   MeanX /= NItem;
   MeanY /= NItem;
   
   /* Calculate correlation coefficient                                 */
   SumXSq = SumYSq = SumXY = (REAL)0.0;
   for(i=0; i<NItem; i++)
   {
      xt = x[i] - MeanX;
      yt = y[i] - MeanY;
      SumXSq += (xt * xt);
      SumYSq += (yt * yt);
      SumXY  += (xt * yt);
   }

   r = SumXY / (REAL)sqrt((double)(SumXSq * SumYSq));
   
   return(r);
}


