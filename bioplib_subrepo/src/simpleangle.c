/*************************************************************************

   Program:    
   File:       simpleangle.c
   
   Version:    V1.5
   Date:       27.03.95
   Function:   
   
   Copyright:  (c) SciTech Software, 1993-5
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

*************************************************************************/
/* Includes
*/
#include <math.h>
#include "MathType.h"
#include "macros.h"

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
/*>REAL simpleangle(REAL ang)
   --------------------------
   Input:   REAL    ang         An angle
   Returns: REAL                Simplified angle
   
   Simplifies a signed angle to an unsigned angle <=2*PI

   07.02.89 Original    By: ACRM
   04.03.91 Fixed return value
   16.06.93 Changed float to REAL
*/
REAL simpleangle(REAL ang)
{
   /* Reduce to less than 360 degrees                                   */
   while(ang > 2*PI) ang -= 2*PI;
   
   if(ang >= 0.0 && ang <= PI)         /* 1st & 2nd quadrant +ve        */
      return(ang);
   else if(ang > PI)                   /* 3rd & 4th quadrant +ve        */
      ang = 2*PI - ang;
   else if(ang < 0.0 && ang > -PI)     /* 1st & 2nd quadrant -ve        */
      ang = ABS(ang);
   else                                /* 3rd & 4th quadrant -ve        */
      ang = 2*PI + ang;
   
   return(ang);
}

