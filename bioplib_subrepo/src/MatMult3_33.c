/*************************************************************************

   Program:    
   File:       MatMult3_33.c
   
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
/*>void MatMult3_33(VEC3F vecin, REAL matin[3][3], VEC3F *vecout)
   -------------------------------------------------------------
   Input:   VEC3F vecin        Vector to be multiplied
            REAL  matin[3][3]  Rotation matrix
   Output:  VEC3F *vecout      Output multiplied vector

   Multiply a 3-vector by a 3x3 matrix

   30.09.92 Original
   03.08.93 Changed multiplication to standard direction
*/
void MatMult3_33(VEC3F vecin, 
                 REAL  matin[3][3], 
                 VEC3F *vecout)
{
   vecout->x = vecin.x * matin[0][0] +
               vecin.y * matin[1][0] +
               vecin.z * matin[2][0];
   vecout->y = vecin.x * matin[0][1] +
               vecin.y * matin[1][1] +
               vecin.z * matin[2][1];
   vecout->z = vecin.x * matin[0][2] +
               vecin.y * matin[1][2] +
               vecin.z * matin[2][2];
}

