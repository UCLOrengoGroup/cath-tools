/***************************************************************************

   Program:    
   File:       cssr.h
   
   Version:    V1.0R
   Date:       09.09.91
   Function:   Defines for CSSR handling
   
   Copyright:  SciTech Software 1991
   Author:     Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      martin@biochem.ucl.ac.uk


****************************************************************************

   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.

   providing it is made clear that this program is free and that the source
   code is provided with the program.

****************************************************************************

   Description:
   ============
   Structure definitions for ReadCSSR()

****************************************************************************

   Usage:
   ======

****************************************************************************

   Revision History:
   =================

***************************************************************************/
#ifndef _CSSR_H
#define _CSSR_H

#include "MathType.h"
#include "pdb.h"

struct cssr_entry
{
   REAL  charge;
   REAL  x,y,z;
   struct cssr_entry *next;
   int   atnum;
   int   group;
   int   link[8];
   char  atnam[8];
} ;

typedef struct cssr_entry CSSR;

#define CLEAR_CSSR(p) p->atnum=0; \
                      strcpy(p->atnam,"    "); \
                      p->group=0; \
                      p->x = 0.0; p->y = 0.0; p->z = 0.0; \
                      p->charge = 0.0; \
                      p->link[0] = p->link[1] = p->link[2] = p->link[3] = 0; \
                      p->link[4] = p->link[5] = p->link[6] = p->link[7] = 0; \
                      p->next = NULL

/************************************************************************/
/* Prototypes
*/
CSSR *ReadCSSR(FILE *fp, int *natom, char *name, char *title);
PDB *ReadCSSRasPDB(FILE *fp, int *natom);
void NormaliseCSSR(CSSR *cssr, REAL cell[3], REAL alpha, REAL beta,
                   REAL gamma);
void NormalisePDB(PDB *pdb, REAL cell[3], REAL alpha, REAL beta,
                  REAL gamma);
void ortho(REAL cell[3], REAL alpha, REAL beta, REAL gamma,
           REAL amatrx[3][3], int isw, int ncode);
void padterm(char *string, int len);
void WriteCSSR(FILE *fp, CSSR *cssr, char *name, char *title);

#endif
