/*************************************************************************

   Program:    
   File:       MatMult33_33.c
   
   Version:    V1.6
   Date:       27.09.95
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
   V1.0  06.09.91 Original
   V1.0a 01.06.92 Documented
   V1.1  30.09.92 Matrix multiplication added
   V1.2  10.06.93 void return from matrix multiplication
   V1.3  22.07.93 Added CreateRotMat()
   V1.4  03.08.93 Changed matrix multiplication to standard direction
   V1.5  28.07.95 Added VecDist()
   V1.6  27.09.95 Added MatMult33_33()

*************************************************************************/
/* Includes
*/
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
/*>void MatMult33_33(REAL a[3][3], REAL b[3][3], REAL out[3][3])
   -------------------------------------------------------------
   Input:   REAL  a[3][3]      Matrix to be multiplied
            REAL  b[3][3]      Matrix to be multiplied
   Output:  REAL  out[3][3]    Output matrix

   Multiply two 3x3 matrices

   27.09.95 Original
*/
void MatMult33_33(REAL a[3][3], REAL b[3][3], REAL out[3][3])
{
   int  i, j, k;
   REAL ab;
   
   for(i=0; i<3; i++)
   {
      for(j=0; j<3; j++)
      {
         ab = (REAL)0.0;
         for(k=0; k<3; k++)
         {
            ab += a[i][k]*b[k][j];
         }
         out[i][j]=ab;
      }
   }
}

         
   
      
   
