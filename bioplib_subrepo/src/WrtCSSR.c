/*************************************************************************

   Program:    
   File:       WrtCSSR.c
   
   Version:    V1.3R
   Date:       01.03.94
   Function:   Write a CSSR file
   
   Copyright:  (c) SciTech Software 1991-4
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
   WriteCSSR(fp,cssr,name,title)
   -----------------------------
   This subroutine will write a CSSR file from a linked list of structures 
   of type cssr_entry. 
   The strucure is set up by including the file "cssr.h". For details of 
   the structure, see this file.

**************************************************************************

   Usage:
   ======
   WriteCSSR(fp,cssr,name,title)
   Input:   FILE     *fp      A pointer to type FILE in which the
                              CSSR file is stored.
            CSSR     *cssr    A pointer to the first allocated item of
                              the CSSR linked list
            char     *name    The molecule's name.
            char     *title   Title on the molecule.

**************************************************************************

   Revision History:
   =================
   V1.0  22.09.91 Original
   V1.1  01.06.92 Autodoc'd. Added FPU check.
   V1.2  10.06.93 void return; float->REAL
   V1.3  27.07.93 %f -> %lf
   V1.3  01.03.94 %lf -> %f  (!)

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "MathType.h"

#include "cssr.h"             /* Note cssr includes pdb                 */
#include "macros.h"

/************************************************************************/
/* Defines and macros
*/
#define CR 13
#define LF 10

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>void WriteCSSR(FILE *fp, CSSR *cssr, char *name, char *title)
   -------------------------------------------------------------
   Input:   FILE     *fp      A pointer to type FILE in which the
                              CSSR file is stored.
            CSSR     *cssr    A pointer to the first allocated item of
                              the CSSR linked list
            char     *name    The molecule's name.
            char     *title   Title on the molecule.

   Write a CSSR file from a CSSR linked list.

   22.09.91 Original
   01.06.92 Autodoc'd
   10.06.93 void return; float->REAL
   27.07.93 %f -> %lf
   01.03.94 %lf -> %f (!)
*/
void WriteCSSR(FILE  *fp,
               CSSR  *cssr,
               char  *name,
               char  *title)
{
   REAL  cell[3],
         alpha, beta, gamma;
   int   natom;
   CSSR  *p;
   
   /* Initilialise cell parameters                                      */
   cell[0] = cell[1] = cell[2] = 1.0;
   alpha = beta = gamma = 90.0;
   
   /*** Record 1, assume no REFCODE                                   ***/
   fprintf(fp,"                                      %8.3f%8.3f%8.3f\n",
           cell[0],cell[1],cell[2]);
   
   /*** Record 2                                                      ***/
   fprintf(fp,"                     %8.3f%8.3f%8.3f\n",
           alpha,beta,gamma);

   /*** Record 3                                                      ***/
   /* Assume orthogonal.
      Count the atoms   
   */
   for(p=cssr, natom=0; p; NEXT(p)) natom++;
   fprintf(fp,"%4d   1 %-60s\n",natom,name);
   
   /*** Record 4                                                      ***/
   /* We'll assume the charges are valid (even if 0.0)                  */
   fprintf(fp,"     2 %-60s\n",title);
   
   /*** Remaining records represent the atoms...                      ***/
   for(p=cssr; p; NEXT(p))
   {
      fprintf(fp,"%4d %-4s  %9.5f %9.5f %9.5f \
%4d%4d%4d%4d%4d%4d%4d%4d %8.4f %2d\n",
              p->atnum,
              p->atnam,
              p->x,
              p->y,
              p->z,
              p->link[0],
              p->link[1],
              p->link[2],
              p->link[3],
              p->link[4],
              p->link[5],
              p->link[6],
              p->link[7],
              p->charge,
              p->group);
   }
}

