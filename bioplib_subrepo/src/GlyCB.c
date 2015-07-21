/*************************************************************************

   Program:    
   File:       GlyCB.c
   
   Version:    V1.0
   Date:       04.01.06
   Function:   Add C-beta atoms to glycines as pseudo-atoms for use
               in orientating residues
   Copyright:  (c) Dr. Andrew C. R. Martin, UCL, 2006
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
   Based on the code from HAddPDB.c adds a C-beta to one or all glycine
   residues.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   04.01.06 V1.0   Original  By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "MathType.h"
#include "SysDefs.h"
#include "pdb.h"
#include "macros.h"

/************************************************************************/
/* Defines and macros
*/
#define BONDLEN 1.524
#define ALPHA   54.7

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>BOOL AddCBtoGly(PDB *pdb)
   -------------------------
   I/O:      PDB   *pdb     The PDB linked list for a Glycine
   Returns:  BOOL           Success?

   Adds a CB atom to a glycine. This is used when one needs to orientate
   a residue in a common frame of reference which makes use of the CB.

   04.01.06 Original   By: ACRM
*/
BOOL AddCBtoGly(PDB *pdb)
{
   PDB *n       = NULL,
       *ca      = NULL,
       *cb      = NULL,
       *c       = NULL,
       *o       = NULL;

   REAL x1,y1,z1,x2,y2,z2,x3,y3,z3,xnew1,ynew1,znew1,
        x21,y21,z21,r21,x23,y23,z23,r23,
        cosa,sina,
        xa,ya,za,xb,yb,zb,xs,ys,zs,
        xab,yab,zab,rab,
        xmin,ymin,zmin,
        xapb,yapb,zapb,rapb,
        xplus,yplus,zplus,
        BondLen = BONDLEN, 
        alpha = ALPHA * PI / 180.0;
/* REAL xnew2,ynew2,znew2;                                              */
   
   if(strncmp(pdb->resnam, "GLY", 3))
      return(FALSE);
   
   if((n  = FindAtomInRes(pdb, "N"))==NULL)
      return(FALSE);
      
   if((ca = FindAtomInRes(pdb, "CA"))==NULL)
      return(FALSE);
   if((c  = FindAtomInRes(pdb, "C"))==NULL)
      return(FALSE);
   if((o  = FindAtomInRes(pdb, "O"))==NULL)
      return(FALSE);

   x1 = n->x;
   y1 = n->y;
   z1 = n->z;
   
   x2 = ca->x;
   y2 = ca->y;
   z2 = ca->z;
   
   x3 = c->x;
   y3 = c->y;
   z3 = c->z;
   
   x21=x2-x1;
   y21=y2-y1;
   z21=z2-z1;
   r21=(REAL)sqrt((double)(x21*x21 + y21*y21 + z21*z21));

   x23=x2-x3;
   y23=y2-y3;
   z23=z2-z3;
   r23=(REAL)sqrt((double)(x23*x23 + y23*y23 + z23*z23));
    
   cosa=(REAL)cos((double)alpha);
   sina=(REAL)sin((double)alpha);

   xa=x21/r21;
   ya=y21/r21;
   za=z21/r21;
   xb=x23/r23;
   yb=y23/r23;
   zb=z23/r23;
   xab=xa-xb;
   yab=ya-yb;
   zab=za-zb;
   rab=(REAL)sqrt((double)(xab*xab+yab*yab+zab*zab));
   xmin=xab/rab;
   ymin=yab/rab;
   zmin=zab/rab;
   xapb=xa+xb;
   yapb=ya+yb;
   zapb=za+zb;
   rapb=(REAL)sqrt((double)(xapb*xapb+yapb*yapb+zapb*zapb));
   xplus=xapb/rapb;
   yplus=yapb/rapb;
   zplus=zapb/rapb;
   xs=yplus*zmin-zplus*ymin;
   ys=zplus*xmin-xplus*zmin;
   zs=xplus*ymin-yplus*xmin;
   xnew1=x2+BondLen*(cosa*xplus-sina*xs);
   ynew1=y2+BondLen*(cosa*yplus-sina*ys);
   znew1=z2+BondLen*(cosa*zplus-sina*zs);

   /* Create a PDB record and initialize it to be the same as the O     */
   if((cb = (PDB *)malloc(sizeof(PDB)))==NULL)
   {
      return(FALSE);
   }
   CopyPDB(cb,o);
   /* Put it into the linked list after the backbone oxygen             */
   cb->next = o->next;
   o->next = cb;
   /* Change it to a CB                                                 */
   strcpy(cb->atnam, "CB  ");
   strcpy(cb->atnam_raw, " CB ");
   /* And set the coordinates                                           */
   cb->x = xnew1;
   cb->y = ynew1;
   cb->z = znew1;

/* This is the position of the hydrogen
   xnew2=x2+BondLen*(cosa*xplus+sina*xs);
   ynew2=y2+BondLen*(cosa*yplus+sina*ys);
   znew2=z2+BondLen*(cosa*zplus+sina*zs);
*/

   return(TRUE);
}


/************************************************************************/
/*>BOOL AddCBtoAllGly(PDB *pdb)
   ----------------------------
   I/O:      PDB   *pdb     The PDB linked list
   Returns:  BOOL           Success?

   Adds a CB atom to all glycines in a PDB linked list. This is used 
   when one needs to orientate a residue in a common frame of reference 
   which makes use of the CB.

   04.01.06 Original   By: ACRM
*/
BOOL AddCBtoAllGly(PDB *pdb)
{
   PDB *p;
   
   for(p=pdb; p!=NULL; p=FindNextResidue(p))
   {
      if(!strncmp(p->resnam, "GLY", 3))
      {
         if(!AddCBtoGly(p))
            return(FALSE);
      }
   }
   return(TRUE);
}


/************************************************************************/
/*>PDB *StripGlyCB(PDB *pdb)
   -------------------------
   I/O:      PDB   *pdb     The PDB linked list
   Returns:  PDB   *        The modified linked list

   Removes all Glycine CB pseudo-atoms added by AddGlyCB()
   The linked list is modified in-place, but the return value
   should be used in case the very first item in the linked list
   is a Gly-CB which will be removed by the code.

   04.01.06 Original   By: ACRM
*/
PDB *StripGlyCB(PDB *pdb)
{
   PDB *p,
       *prev = NULL;
   
   for(p=pdb; p!=NULL; )
   {
      if(!strncmp(p->resnam, "GLY", 3))
      {
         if(!strncmp(p->atnam, "CB  ", 4))
         {
            if(prev!=NULL)
            {
               prev->next = p->next;
               free(p);
               p=prev->next;
            }
            else
            {
               PDB *q;
               q=p;
               NEXT(p);
               free(q);
               pdb = p;
            }
         }
         else
         {
            prev = p;
            NEXT(p);
         }
      }
      else
      {
         prev=p;
         NEXT(p);
      }
   }
   return(pdb);
}
