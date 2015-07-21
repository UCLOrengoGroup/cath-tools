/*************************************************************************

   Program:    
   File:       DistPtVect.c
   
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
#include "MathUtil.h"

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
/*>REAL DistPtVect(VEC3F Point, VEC3F End1, VEC3F End2)
   ----------------------------------------------------
   Input:   VEC3F   Point     The coordinates of a point
            VEC3F   End1      Coordinates of one end of vector
            VEC3F   End2      Coordinates of other end of vector
   Returns: REAL              The distance from pt to line

   Calculate the distance from a point to a vector described by two
   end points

   18.06.96 Original   By: ACRM
*/
REAL DistPtVect(VEC3F Point, VEC3F End1, VEC3F End2)
{
   VEC3F Vec,
         UVec,
         PQVec,
         PRVec;
   REAL  len;
   
   /* Find the vector from End1 to End2                                 */
   VecSub3(&Vec, End2, End1);

   /* Find the length of this vector                                    */
   len = VecLen3(Vec);

   /* Now calculate the unit vector                                     */
   UVec.x = Vec.x / len;
   UVec.y = Vec.y / len;
   UVec.z = Vec.z / len;

   /* Calculate the vector from the point to an arbitrary point on the
      line (we'll choose End1)
   */
   PQVec.x = End1.x - Point.x;
   PQVec.y = End1.y - Point.y;
   PQVec.z = End1.z - Point.z;

   /* PRVect is the cross product of PQVect with the unit vector        */
   CrossProd3(&PRVec, PQVec, UVec);

   /* The length we want is the length of this vector                   */
   return(VecLen3(PRVec));
}
