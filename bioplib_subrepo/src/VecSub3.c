/*************************************************************************

   Program:    
   File:       VecSub3.c
   
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
/*>void VecSub3(VEC3F *Out, VEC3F In1, VEC3F In2)
   ----------------------------------------------
   Input:   VEC3F   In1       First vector
            VEC3F   In2       Second vector
   Output:  VEC3F   Out       Output vector

   Subtract 2 vectors
   18.06.96 Original   By: ACRM
*/
void VecSub3(VEC3F *Out, VEC3F In1, VEC3F In2)
{
   Out->x = In1.x - In2.x;
   Out->y = In1.y - In2.y;
   Out->z = In1.z - In2.z;
}


