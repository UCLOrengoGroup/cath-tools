/*************************************************************************

   Program:    
   File:       SetChi.c
   
   Version:    V1.2
   Date:       01.03.94
   Function:   
   
   Copyright:  (c) Dr. Andrew C. R. Martin, University of Reading, 2002
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
   V1.1  01.03.94
   V1.2  27.02.98 Removed unreachable break from switch()

*************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include <math.h>

#include "MathType.h"
#include "pdb.h"
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
/*>void SetChi(PDB *pdb, PDB *next, REAL chi, int type)
   ----------------------------------------------------
   I/O:    PDB  *pdb    PDB linked list to change
   Input:  PDB  *next   If NULL, move all atoms in the linked list from
                        the last atom in the torsion. Otherwise move
                        atoms up to (but not including) next. Normally
                        this would be the start of the next residue.
           REAL chi     Sidechain chi angle to set
           int  type    Torsion angle to set (as defined below)

   Sets a sidechain torsion angle in a pdb linked list. The routine 
   assumes standard atom ordering: N,CA,C,O,s/c with standard order in
   the s/c.
   The type input parameter is defined as follows:
   
         type     Atom names        Sequential atom numbers
         --------------------------------------------------
         0        N,  CA, CB, XG    (0 - 1 - 4 - 5)
         1        CA, CB, XG, XD    (1 - 4 - 5 - 6)
         2        CB, XG, XD, XE    (4 - 5 - 6 - 7)
         3        XG, XD, XE, XZ    (5 - 6 - 7 - 8)
   
   13.05.92 Original
   27.02.98 Removed unreachable break from switch()
*/
void SetChi(PDB   *pdb,
            PDB   *next, 
            REAL  chi, 
            int   type)
{
   int   natoms,
         nmove,
         j;
   REAL  *x = NULL,
         *y = NULL,
         *z = NULL,
         s,
         n1, n2, n3,
         sinrot,
         cosrot,
         matrix[3][3],
         CurrentChi;
   PDB   *one,
         *two,
         *three,
         *four,
         *p;
   VEC3F base;
   
   /* First count number of atoms                                       */
   for(p=pdb, natoms = 0; p!=next; NEXT(p), natoms++) ;
   
   /* Allocate memory for the coordinate lists                          */
   x = (REAL *)malloc(natoms * sizeof(REAL));
   y = (REAL *)malloc(natoms * sizeof(REAL));
   z = (REAL *)malloc(natoms * sizeof(REAL));
   if(x==NULL || y==NULL || z==NULL) goto Cleanup;
   
   /* Find the current value of the torsion angle                       */
   CurrentChi = CalcChi(pdb, type);
   
   /* Calc rotation required.                                           */
   chi -= CurrentChi;
   
   /* Get pointers to the appropriate atoms                             */
   switch(type)
   {
   case 0:              /* N,  CA, CB, XG    (0 - 1 - 4 - 5)            */
      one   = GetPDBByN(pdb, 0);
      two   = GetPDBByN(pdb, 1);
      three = GetPDBByN(pdb, 4);
      four  = GetPDBByN(pdb, 5);
      break;
   case 1:              /* CA, CB, XG, XD    (1 - 4 - 5 - 6)            */
      one   = GetPDBByN(pdb, 1);
      two   = GetPDBByN(pdb, 4);
      three = GetPDBByN(pdb, 5);
      four  = GetPDBByN(pdb, 6);
      break;
   case 2:              /* CB, XG, XD, XE    (4 - 5 - 6 - 7)            */
      one   = GetPDBByN(pdb, 4);
      two   = GetPDBByN(pdb, 5);
      three = GetPDBByN(pdb, 6);
      four  = GetPDBByN(pdb, 7);
      break;
   case 3:              /* XG, XD, XE, XZ    (5 - 6 - 7 - 8)            */
      one   = GetPDBByN(pdb, 5);
      two   = GetPDBByN(pdb, 6);
      three = GetPDBByN(pdb, 7);
      four  = GetPDBByN(pdb, 8);
      break;
   default:
      return;
   }
   
   /* Copy the mobile atoms into the arrays                             */
   x[0] = two->x;
   y[0] = two->y;
   z[0] = two->z;
   x[1] = three->x;
   y[1] = three->y;
   z[1] = three->z;
   for(p=four, nmove=2; p!=next && nmove<natoms; NEXT(p), nmove++)
   {
      x[nmove] = p->x;
      y[nmove] = p->y;
      z[nmove] = p->z;
   }
   
   base.x = x[0];
   base.y = y[0];
   base.z = z[0];
   
   /* Move the atoms from atom two on to the origin                     */
   for(j=0; j<nmove; j++)
   {
      x[j] -= base.x;
      y[j] -= base.y;
      z[j] -= base.z;
   }
   
   /* RotatePDB about an arbitrary axis using routine C-11 from Rogers & 
      Adams, `Mathematical Elements for Computer Graphics'
   */
   s = (REAL)sqrt((double)(x[1] * x[1] + y[1] * y[1] + z[1] * z[1]));
   
   n1 = x[1]/s;
   n2 = y[1]/s;
   n3 = z[1]/s;
   
   cosrot = (REAL)cos((double)chi);
   sinrot = (REAL)sin((double)chi);
   
   /* Set up the transformation matrix                                  */
   matrix[0][0] = n1*n1+(1-n1*n1)*cosrot;
   matrix[0][1] = n1*n2*(1-cosrot)+n3*sinrot;
   matrix[0][2] = n1*n3*(1-cosrot)-n2*sinrot;
   matrix[1][0] = n1*n2*(1-cosrot)-n3*sinrot;
   matrix[1][1] = n2*n2+(1-n2*n2)*cosrot;
   matrix[1][2] = n2*n3*(1-cosrot)+n1*sinrot;
   matrix[2][0] = n1*n3*(1-cosrot)+n2*sinrot;
   matrix[2][1] = n2*n3*(1-cosrot)-n1*sinrot;
   matrix[2][2] = n3*n3+(1-n3*n3)*cosrot;
   
   /* Do the matrix multiplication                                      */
   for(j=0; j<nmove; j++)
   {
      REAL  tempx,
            tempy,
            tempz;
            
      tempx =  x[j] * matrix[0][0] +
               y[j] * matrix[1][0] +
               z[j] * matrix[2][0];
      tempy =  x[j] * matrix[0][1] +
               y[j] * matrix[1][1] +
               z[j] * matrix[2][1];
      tempz =  x[j] * matrix[0][2] +
               y[j] * matrix[1][2] +
               z[j] * matrix[2][2];

      x[j] = tempx;
      y[j] = tempy;
      z[j] = tempz;
   }

   /* TranslatePDB back from the origin                                 */
   for(j=0; j<nmove; j++)
   {
      x[j] += base.x;
      y[j] += base.y;
      z[j] += base.z;
   }
   
   /* Copy the new coordinates back                                     */
   two->x   = x[0];
   two->y   = y[0];
   two->z   = z[0];
   three->x = x[1];
   three->y = y[1];
   three->z = z[1];
   for(p=four, nmove=2; p!=next && nmove<natoms; NEXT(p), nmove++)
   {
      p->x = x[nmove];
      p->y = y[nmove];
      p->z = z[nmove];
   }
   
Cleanup:
   if(x) free(x);
   if(y) free(y);
   if(z) free(z);
}


