/*************************************************************************

   Program:    
   File:       CreateRotMat.c
   
   Version:    V1.6R
   Date:       27.09.95
   Function:   Simple matrix and vector operations
   
   Copyright:  (c) SciTech Software 1991-5
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      martin@biochem.ucl.ac.uk
               
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
/*>void CreateRotMat(char direction, REAL angle, REAL matrix[3][3])
   ----------------------------------------------------------------
   Input:   char direction    Axis about which to rotate
            REAL angle        Angle (in rads) to rotate
   Output:  REAL matrix[3][3] Rotation matrix

   Create a 3x3 rotation matrix. Takes a direction as a single character
   ('x', 'y', or 'z'), an angle (in rads) and outputs a rotation matrix
   
   22.07.93 Original    By: ACRM
*/
void CreateRotMat(char direction, REAL angle, REAL matrix[3][3])
{
   int   i, j,
         m, 
         m1,
         m2;
   REAL  CosTheta,
         SinTheta;
         
   /* Initialise matrix to all 0.0                                      */
   for(i=0; i<3; i++)
      for(j=0; j<3; j++)
         matrix[i][j] = 0.0;
   
   /* Select the items that need to be filled in                        */
   switch(direction)
   {
   case 'x':   case 'X':
      m = 0;
      break;
   case 'y':   case 'Y':
      m = 1;
      break;
   case 'z':   case 'Z':
      m = 2;
      break;
   default:                   /* Just return the unit matrix            */
      for(i=0; i<3; i++)
         matrix[i][i] = 1.0;
      return;
   }
   
   /* Find which items these relate to                                  */
   m1 = (m+1)  % 3;
   m2 = (m1+1) % 3;
   
   /* Fill in the values                                                */
   matrix[m][m]   = 1.0;
   CosTheta       = (REAL)cos((double)angle);
   SinTheta       = (REAL)sin((double)angle);
   matrix[m1][m1] = CosTheta;
   matrix[m2][m2] = CosTheta;
   matrix[m1][m2] = SinTheta;
   matrix[m2][m1] = -SinTheta;
}


