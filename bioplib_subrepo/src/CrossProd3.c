/*************************************************************************

   Program:    
   File:       CrossProd3.c
   
   Version:    V1.2
   Date:       06.10.98
   Function:   General maths/stats/vector functions
   
   Copyright:  (c) SciTech Software 1996-8
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
/*>void CrossProd3(VEC3F *Out, VEC3F In1, VEC3F In2)
   -------------------------------------------------
   Input:   VEC3F   In1       First vector
            VEC3F   In2       Second vector
   Output:  VEC3F   Out       Output vector

   Calculate the cross product of 2 vectors
   18.06.96 Original   By: ACRM
*/
void CrossProd3(VEC3F *Out, VEC3F In1, VEC3F In2)
{
   Out->x = In1.y*In2.z - In1.z*In2.y;
   Out->y = In1.z*In2.x - In1.x*In2.z;
   Out->z = In1.x*In2.y - In1.y*In2.x;
}

