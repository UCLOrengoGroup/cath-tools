/*************************************************************************

   Program:    
   File:       TrueAngle.c
   
   Version:    V1.5
   Date:       27.03.95
   Function:   
   
   Copyright:  (c) SciTech Software 1993-5
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
   REAL TrueAngle(REAL opp, REAL adj)
   Input:   REAL     opp         Length of opposite side
            REAL     adj         Length of adjacent side
   Returns: REAL                 The angle from 0 to 2PI

   Returns the true positive angle between 0 and 2PI given the opp and
   adj lengths

**************************************************************************

   Revision History:
   =================
   V1.0  07.02.91 Original
   V1.1  17.02.91 Corrected comments to new standard and added phi()
   V1.2  04.03.91 angle() and phi() now return _correct_ values!
   V1.3  01.06.92 ANSIed
   V1.4  08.12.92 Changed abs() to ABS() from macros.h
   V1.5  27.03.95 Added TrueAngle()

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
/*>REAL TrueAngle(REAL opp, REAL adj)
   ----------------------------------
   Input:   REAL     opp     Opposite length
            REAL     adj     Adjacent length
   Returns: REAL             Angle between 0 and 2PI

   Return the +ve angle between 0 and 2PI given the opp and adj values.

   25.07.94 Original    By: ACRM
*/
REAL TrueAngle(REAL opp, REAL adj)
{
   REAL ang;
   
   if(adj != 0.0)
   {
      ang = (REAL)atan((double)(opp/adj));

      /* 4th quadrant; ang -ve so add 2PI                             */
      if(opp < 0.0 && adj > 0.0) ang += 2*PI;

      /* 2nd & 3rd quadrant; add PI                                     */
      if(adj < 0.0) ang += PI;
   }
   else
   {
      if(opp>0.0)                /* 1st->2nd quadrant boundary          */
         ang = PI/2.0;
      else                       /* 3rd->4th quadrant boundary          */
         ang = 3.0*PI/2.0;
   }
   
   if(opp == 0.0)
   {
      if(adj > 0.0)              /* 4th->1st quadrant boundary          */
         ang = 0.0;
      else                       /* 2nd->3rd quadrant boundary          */
         ang = PI;
   }

   return(ang);
}

