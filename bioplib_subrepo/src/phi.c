/*************************************************************************

   Program:    
   File:       phi.c
   
   Version:    V1.5
   Date:       27.03.95
   Function:   Calculate a torsion angle
   
   Copyright:  (c) SciTech Software 1993-5
   Author:     Dr. Andrew C. R. Martin
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
/*>REAL phi(REAL xi,REAL yi,REAL zi,REAL xj,REAL yj,REAL zj,
            REAL xk,REAL yk,REAL zk,REAL xl,REAL yl,REAL zl)
   ---------------------------------------------------------
   Input:   REAL    xi,yi,zi    Input coordinates
                    xj,yj,zj
                    xk,yk,zk
                    xl,yl,zl
   Returns: REAL                The torsion angle between the 4 atoms

   Calculates the torsion angle described by 4 sets of coordinates.

   04.03.91 Original    By: ACRM
   16.06.93 Changed float to REAL
*/
REAL phi(REAL xi,
         REAL yi,
         REAL zi,
         REAL xj,
         REAL yj,
         REAL zj,
         REAL xk,
         REAL yk,
         REAL zk,
         REAL xl,
         REAL yl,
         REAL zl)
{
   REAL xij,yij,zij,
        xkj,ykj,zkj,
        xkl,ykl,zkl,
        dxi,dyi,dzi,
        gxi,gyi,gzi,
        bi,bk,ct,
        boi2,boj2,
        z1,z2,ap,s,
        bioj,bjoi;


   /* Calculate the vectors C,B,C                                       */
   xij = xi - xj;
   yij = yi - yj;
   zij = zi - zj;
   xkj = xk - xj;
   ykj = yk - yj;
   zkj = zk - zj;
   xkl = xk - xl;
   ykl = yk - yl;
   zkl = zk - zl;

   /* Calculate the normals to the two planes n1 and n2
      this is given as the cross products:
       AB x BC
      --------- = n1
      |AB x BC|

       BC x CD
      --------- = n2
      |BC x CD|
   */
   dxi = yij * zkj - zij * ykj;     /* Normal to plane 1                */
   dyi = zij * xkj - xij * zkj;
   dzi = xij * ykj - yij * xkj;
   gxi = zkj * ykl - ykj * zkl;     /* Normal to plane 2                */
   gyi = xkj * zkl - zkj * xkl;
   gzi = ykj * xkl - xkj * ykl;

   /* Calculate the length of the two normals                           */
   bi = dxi * dxi + dyi * dyi + dzi * dzi;
   bk = gxi * gxi + gyi * gyi + gzi * gzi;
   ct = dxi * gxi + dyi * gyi + dzi * gzi;

   boi2 = 1./bi;
   boj2 = 1./bk;
   bi   = (REAL)sqrt((double)bi);
   bk   = (REAL)sqrt((double)bk);

   z1   = 1./bi;
   z2   = 1./bk;
   bioj = bi * z2;
   bjoi = bk * z1;
   ct   = ct * z1 * z2;
   if (ct >  1.0)   ct = 1.0;
   if (ct < (-1.0)) ct = -1.0;
   ap   = acos(ct);

   s = xkj * (dzi * gyi - dyi * gzi)
     + ykj * (dxi * gzi - dzi * gxi)
     + zkj * (dyi * gxi - dxi * gyi);

   if (s < 0.0) ap = -ap;

   ap = (ap > 0.0) ? PI-ap : -(PI+ap);

   return(ap);
}

