/*************************************************************************

   Program:    
   File:       pldist.c
   
   Version:    V1.0
   Date:       16.11.99
   Function:   Calculate distance from a point to a line
   
   Copyright:  (c) University of Reading / Dr. Andrew C. R. Martin 1999
   Author:     Dr. Andrew C. R. Martin
   Address:    School of Animal and Microbial Sciences,
               The University of Reading,
               Whiteknights,
               P.O. Box 228,
               Reading RG6 6AJ.
               England.
   Phone:      +44 (0)118 987 5123 Extn. 7022
   Fax:        +44 (0)118 931 0180
   EMail:      a.c.r.martin@reading.ac.uk
               andrew@stagleys.demon.co.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

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
#include <stdio.h>
#include <math.h>
#include "MathType.h"

/************************************************************************/
/* Defines and macros
*/
#define CROSSPRODUCT(p1,p2,p3) \
        (p3).x = (p1).y*(p2).z - (p1).z*(p2).y; \
        (p3).y = (p1).z*(p2).x - (p1).x*(p2).z; \
        (p3).z = (p1).x*(p2).y - (p1).y*(p2).x
#define DOTPRODUCT(v1,v2) ((v1).x*(v2).x + (v1).y*(v2).y + (v1).z*(v2).z)

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>REAL PointLineDistance(REAL Px, REAL Py, REAL Pz,
                          REAL P1x, REAL P1y, REAL P1z,
                          REAL P2x, REAL P2y, REAL P2z,
                          REAL *Rx, REAL *Ry, REAL *Rz,
                          REAL *frac)
   ----------------------------------------------------
   Input:    REAL   Px          Point x coordinate
             REAL   Py          Point y coordinate
             REAL   Pz          Point z coordinate
             REAL   P1x         Line start x coordinate
             REAL   P1y         Line start y coordinate
             REAL   P1z         Line start z coordinate
             REAL   P2x         Line end x coordinate
             REAL   P2y         Line end y coordinate
             REAL   P2z         Line end z coordinate
    Output:  REAL   *Rx         Nearest point on line x coordinate
             REAL   *Ry         Nearest point on line y coordinate
             REAL   *Rz         Nearest point on line z coordinate
             REAL   *frac       Fraction along P1-P2 of R
    Returns: REAL               Distance from P to R

    Calculates the shortest distance from a point P to a line between
    points P1 and P2. This value is returned.

    If the Rx,Ry,Rz pointers are all non-NULL, then the point on the
    line nearest to P is output.

    If the frac pointer is non-NULL then the fraction of R along the
    P1-P2 vector is output. Thus:
          R==P1 ==> frac=0
          R==P2 ==> frac=1
    Thus if (0<=frac<=1) then the point R is within the line segment
    P1-P2

    16.11.99 Original   By: ACRM
*/
REAL PointLineDistance(REAL Px, REAL Py, REAL Pz,
                       REAL P1x, REAL P1y, REAL P1z,
                       REAL P2x, REAL P2y, REAL P2z,
                       REAL *Rx, REAL *Ry, REAL *Rz,
                       REAL *frac)
{
   VEC3F A, u, Q, PQ, PR, QP, QP2;
   REAL  alen, len, f;
   
   
   /* Calculate vector from P1 to P2                                    */
   A.x = P2x - P1x;
   A.y = P2y - P1y;
   A.z = P2z - P1z;
   
   /* Calculate length of this vector                                   */
   alen = sqrt(DOTPRODUCT(A,A));

   /* If the two ends of the line are coincident then return the distance
      from either of them
   */
   if(alen==(REAL)0.0)
   {
      len = sqrt((Px-P1x)*(Px-P1x) + 
                 (Py-P1y)*(Py-P1y) + 
                 (Pz-P1z)*(Pz-P1z));
      if(frac!=NULL)
         *frac = 0.0;
      if(Rx != NULL && Ry != NULL && Rz != NULL)
      {
         *Rx = P1x;
         *Ry = P1y;
         *Rz = P1z;
      }

      return(len);
   }

   /* Calculate the unit vector along A                                 */
   u.x = A.x / alen;
   u.y = A.y / alen;
   u.z = A.z / alen;
   
   /* Select Q as any point on A, we'll make it P1                      */
   Q.x = P1x;
   Q.y = P1y;
   Q.z = P1z;
   
   /* Calculate vector PQ                                               */
   PQ.x = Q.x - Px;
   PQ.y = Q.y - Py;
   PQ.z = Q.z - Pz;

   /* Vector PR is the cross product of PQ and the unit vector
      along A (i.e. u)
   */
   CROSSPRODUCT(PQ, u, PR);
   
   /* And the length of that vector is the length we want               */
   len = sqrt(DOTPRODUCT(PR,PR));


   if(frac != NULL || Rx != NULL || Ry != NULL || Rz != NULL)
   {
      /*** OK we now know how far the point is from the line, so we   ***
       *** now want to calculate where the closest point (R) on the   ***
       *** line is to point P                                         ***/

      /* Find the projection of QP onto QP2                             */
      QP.x = Px - Q.x;
      QP.y = Py - Q.y;
      QP.z = Pz - Q.z;
      
      QP2.x = P2x - Q.x;
      QP2.y = P2y - Q.y;
      QP2.z = P2z - Q.z;
      
      f = DOTPRODUCT(QP, QP2) / sqrt(DOTPRODUCT(QP2, QP2));
      if(frac != NULL)
      {
         *frac = f/alen;
      }
      
      /* Find point R: this is the fraction f of the unit vector along 
         P1-P2 added onto Q 
      */
      if(Rx != NULL && Ry != NULL && Rz != NULL)
      {
         *Rx = Q.x + f * u.x;
         *Ry = Q.y + f * u.y;
         *Rz = Q.z + f * u.z;
      }
   }
   
   return(len);
}


/************************************************************************/
#ifdef DEMO
int main(int argc, char **argv)
{
   REAL Px, Py, Pz,
        P1x, P1y, P1z,
        P2x, P2y, P2z,
        Rx, Ry, Rz, d, f;
 
   Px = 5;   Py = 2;   Pz = 0;
   P1x = 5;  P1y = 2;  P1z = 0;
   P2x = 10; P2y = 2;  P2z = 0;
   
   d = PointLineDistance(Px, Py, Pz,
                         P1x, P1y, P1z,
                         P2x, P2y, P2z,
                         &Rx, &Ry, &Rz, &f);
   
   printf("*** Distance is %f; Point R is %f %f %f; f is %f\n",
          d,Rx,Ry,Rz,f);
   
   return(0);
}
#endif
