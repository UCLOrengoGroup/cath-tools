/*************************************************************************

   Program:    
   File:       CalcCellTrans.c
   
   Version:    V1.0R
   Date:       12.10.95
   Function:   Calculate offsets for creating a crystal lattice
   
   Copyright:  (c) SciTech Software 1993-1995
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
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
   V1.0R 12.10.05 Original

*************************************************************************/
/* Defines required for includes
*/

/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>

#include "SysDefs.h"
#include "MathType.h"
#include "pdb.h"
#include "fsscanf.h"

/************************************************************************/
/* Defines
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>void CalcCellTrans(VEC3F UnitCell, VEC3F CellAngles, 
                      VEC3F *xtrans, VEC3F *ytrans, VEC3F *ztrans)
   ---------------------------------------------------------------
   Input:   VEC3F  UnitCell       The unit cell dimensions
            VEC3F  CellAngles     The unit cell angles
   Output:  VEC3F  *xtrans        Translation to apply along X axis
            VEC3F  *ytrans        Translation to apply along Y axis
            VEC3F  *ztrans        Translation to apply along Z axis

   Calculates the offsets to apply in X, Y and Z directions for creating
   a crystal lattice from the unit cell parameters.

   11.10.95 Original    By: ACRM, Based on code from Rasmol by Roger
                        Sayle (ros@dcs.ed.ac.uk, ras32425@ggr.co.uk)
*/
void CalcCellTrans(VEC3F UnitCell, VEC3F CellAngles, 
                   VEC3F *xtrans, VEC3F *ytrans, VEC3F *ztrans)
{
   REAL lena, lenb, lenc,
        cosa, cosb, cosg, sing,
        tmpx, tmpy, tmpz, temp;
   
   lena =  UnitCell.x;                         /* A                     */
   lenb =  UnitCell.y;                         /* B                     */
   lenc = -1.0*UnitCell.z;                     /* C                     */

   cosa = (REAL)cos((double)(CellAngles.x));   /* Alpha                 */
   cosb = (REAL)cos((double)(CellAngles.y));   /* Beta                  */
   cosg = (REAL)cos((double)(CellAngles.z));   /* Gamma                 */
   sing = (REAL)sin((double)(CellAngles.z));   /* Gamma                 */

   temp = cosa*cosa + cosb*cosb + cosg*cosg - 2.0*cosa*cosb*cosg;
   tmpx = cosb; 
   tmpy = (cosa - cosb*cosg)/sing;
   tmpz = (REAL)sqrt((double)(1.0-temp))/sing;
   
   xtrans->x = lena;
   xtrans->y = (REAL)0.0;
   xtrans->z = (REAL)0.0;
   
   ytrans->x = lenb*cosg;
   ytrans->y = lenb*sing;
   ytrans->z = (REAL)0.0;
   
   ztrans->x = lenc*tmpx;
   ztrans->y = lenc*tmpy;
   ztrans->z = lenc*tmpz;
}


