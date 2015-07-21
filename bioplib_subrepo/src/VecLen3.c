/*************************************************************************

   Program:    
   File:       VecLen3.c
   
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
/*>REAL VecLen3(VEC3F Vec)
   -----------------------
   Input:   VEC3F   Vec       Vector
   Returns: REAL              Length of vector

   18.06.96 Original   By: ACRM
*/
REAL VecLen3(VEC3F Vec)
{
   return((REAL)sqrt((double)((Vec.x * Vec.x) + 
                              (Vec.y * Vec.y) + 
                              (Vec.z * Vec.z))));
}


