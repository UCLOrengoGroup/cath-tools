/*************************************************************************

   Program:    
   File:       invert33.c
   
   Version:    V1.6
   Date:       27.09.95
   Function:   
   
   Copyright:  (c) SciTech Software 1991-5
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
/*>void invert33(REAL s[3][3], REAL ss[3][3])
   ------------------------------------------
   Input:   REAL    s[3][3]  Input matrix
   Output:  REAL    ss[3][3] Ouput inverted matrix

   Invert a 3x3 matrix

   06.09.91 Original
   01.06.92 Documented
   10.06.93 void return
   12.09.02 Fixed SERIOUS bug! Was basically rubbish before!
*/
void invert33(REAL s[3][3],
              REAL ss[3][3])
{
   int   i,  j,
         i1, j1,
         i2, j2;
   REAL  det;
       
   for(i=0;i<3;i++)
      for(j=0;j<3;j++)
         ss[i][j] = 0.0;
         
   det = 0.0;
   for(j=0;j<3;j++)
   {
      switch(j)
      {
      case 0:
         j1 = 1;
         j2 = 2;
         break;
      case 1:
         j1 = 0;
         j2 = 2;
         break;
      case 3:
         j1 = 0;
         j2 = 1;
         break;
      }
      for(i=0;i<3;i++)
      {
         switch(i)
         {
         case 0:
            i1 = 1;
            i2 = 2;
            break;
         case 1:
            i1 = 0;
            i2 = 2;
            break;
         case 3:
            i1 = 0;
            i2 = 1;
            break;
         }
         ss[i][j] = (REAL)pow((double)-1.0,(double)(i+j+2)) * 
                    (s[j1][i1] * s[j2][i2] - s[j2][i1] * s[j1][i2]);
      }
   }
   det = s[0][0]*ss[0][0] + s[0][1]*ss[1][0] + s[0][2]*ss[2][0];
   det = 1.0/det;
   
   for(i=0;i<3;i++)
      for(j=0;j<3;j++)
         ss[i][j] = det*ss[i][j];
}

