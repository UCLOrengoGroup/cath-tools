/*************************************************************************

   Program:    
   File:       ReadCSSR.c
   
   Version:    V1.5R
   Date:       30.05.02
   Function:   Read a CSSR file
   
   Copyright:  (c) SciTech Software 1991-2002
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
   ReadCSSR(fp,cssr,natom,name,title)
   ----------------------------------
   This subroutine will read a CSSR file of any size and form a linked list
   of the structure.
   This list is contained in a linked set of structures of type cssr_entry. 
   The strucure is set up by including the file "cssr.h". For details of 
   the structure, see this file.

   To define a structure list in which to store the protein, the
   user need only include the file "cssr.h", declare a pointer to a
   structure of type CSSR using the statement:
       CSSR *mycssr;


   ReadCSSRasPDB(fp,cssr,natom)
   ----------------------------
   As ReadCSSR(), but reads the structure into a PDB type linked list.
   Atom connection information is ignored, and charges are placed into the
   B-val column.

   NormaliseCSSR(cssr,cell,alpha,beta,gamma)
   -----------------------------------------
   Given the unit cell dimensions, converts CSSR to orthonormal 
   coordinates.

   NormalisePDB(pdb,cell,alpha,beta,gamma)
   ---------------------------------------
   Given the unit cell dimensions, converts PDB to orthonormal 
   coordinates.

   ortho(cell,alpha,beta,gamma,amatrx,isw,ncode)
   ---------------------------------------------
   Calculates the 3x3 matrix required to convert between fractional and
   orthonormal coordinates given the unit cell dimensions.

**************************************************************************

   Usage:
   ======
   ReadCSSR(fp,cssr,natom,name,title)
   Input:   FILE     *fp      A pointer to type FILE in which the
                              CSSR file is stored.
            CSSR     *cssr    A pointer to the first allocated item of
                              the CSSR linked list
   Output:  int      *natom   Number of atoms read.
            char     *name    The molecule's name.
            char     *title   Title on the molecule.

   ReadCSSRasPDB(fp,pdb,natom)
                              into a PDB linked list
   Input:   FILE     *fp      A pointer to type FILE in which the
                              CSSR file is stored.
            PDB      *pdb     A pointer to the first allocated item of
                              the PDB linked list
   Output:  int      *natom   Number of atoms read.

   NormaliseCSSR(cssr,cell,alpha,beta,gamma)
   I/O:     CSSR     *cssr    Pointer to CSSR linked list
   Input:   REAL     cell[3]  Unit cell dimensions
                     alpha    Unit cell angles
                     beta
                     gamma

   NormalisePDB(pdb,cell,alpha,beta,gamma)
   I/O:     PDB      *pdb     Pointer to PDB linked list
   Input:   REAL     cell[3]  Unit cell dimensions
                     alpha    Unit cell angles
                     beta
                     gamma

   ortho(cell,alpha,beta,gamma,amatrx,isw,ncode)
   Input:   REAL     cell[3]  Unit cell dimensions
                     alpha    Unit cell angles
                     beta
                     gamma
   Output:  REAL     amatrx[3][3]   Returned conversion matrix
   Input:   int      isw      0: Frac-->Ortho,  1: Ortho-->Frac
            int      ncode    Orientation of reciprocal axes wrt true 
                              axes.
                              1: a || xo,          c* || zo
                              2: b || xo,          a* || zo
                              3: c || xo,          b* || zo
                              4: hex a & b || xo,  c* || zo
                              5: a* || xo,         c  || zo

**************************************************************************

   Revision History:
   =================
   V1.0  06.09.91 Original
   V1.1  24.01.92 Fixed for reading files with blank link numbers.
   V1.2  01.06.92 Documented and added FPU check
   V1.3  10.06.93 Tidied for book
   V1.4  27.07.93 Changed I/O to double precision
   V1.5  30.05.02 Changed PDB field from 'junk' to 'record_type'

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "SysDefs.h"
#include "MathType.h"

#include "cssr.h"        /* This includes pdb.h                         */
#include "matrix.h"
#include "general.h"
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
/*>CSSR *ReadCSSR(FILE *fp, int *natom, char *name, char *title)
   -------------------------------------------------------------
   Input:   FILE     *fp      A pointer to type FILE in which the
                              CSSR file is stored.
            CSSR     *cssr    A pointer to the first allocated item of
                              the CSSR linked list
   Output:  int      *natom   Number of atoms read.
            char     *name    The molecule's name.
            char     *title   Title on the molecule.
   Read a CSSR file into a CSSR linked list

   06.09.91 Original
   24.01.92 Fixed for blank link columns (V1.1)
   01.06.92 Documented
   10.06.93 Returns TRUE or FALSE to indicate success
   09.07.93 Changed allocation scheme. Now returns pointer to start of
            list. No need to call init_cssr
   13.07.93 Returns NULL if allocation failed
   27.07.93 Changed I/O to double precision
*/
CSSR *ReadCSSR(FILE  *fp,
               int   *natom,
               char  *name,
               char  *title)
{
   char  buffer[160],
         code[8];
   REAL  cell[3],
         alpha, beta, gamma;
   int   orthonormal,
         nocharges,
         i;
   CSSR  *cssr = NULL,
         *p;
   
   /* Initilialise cell parameters, etc.                                */
   cell[0]  = cell[1]  = cell[2]  = 1.0;
   alpha    = beta     = gamma    = 90.0;
   name[0]  = title[0] = '\0';
   
   *natom   = 0;
   
   /*** Record 1                                                      ***/
   if(!fgets(buffer,160,fp)) return(NULL);
   TERMINATE(buffer);
   
   if(strlen(buffer) >= 28)
   {
      /* We can get a code name                                         */
      strncpy(code,buffer+22,6);
      code[6] = '\0';
   }
   
   if(strlen(buffer) >= 62)
   {
      /* We can get unit cell parameters                                */
      sscanf(buffer+38,"%lf%lf%lf",&cell[0],&cell[1],&cell[2]);
   }
   
   /*** Record 2                                                      ***/
   if(!fgets(buffer,160,fp)) return(NULL);
   TERMINATE(buffer);
   
   if(strlen(buffer) >= 45)
   {
      /* We can get cell angles                                         */
      sscanf(buffer+21,"%lf%lf%lf",&alpha,&beta,&gamma);
   }
   
   /*** Record 3                                                      ***/
   if(!fgets(buffer,160,fp)) return(NULL);
   TERMINATE(buffer);
   
   sscanf(buffer,"%d%d",natom,&orthonormal);
   strcpy(name,buffer+9);
   
   /*** Record 4                                                      ***/
   if(!fgets(buffer,160,fp)) return(NULL);
   TERMINATE(buffer);
   
   sscanf(buffer+4,"%d",&nocharges);
   if(nocharges != 1) nocharges = 0;
   strcpy(title,buffer+7);
   
   /*** Remaining records represent the atoms...                      ***/
   for(i=0; i < *natom; i++)
   {
      /* Read atom record                                               */
      if(!fgets(buffer,160,fp)) return(cssr);
      TERMINATE(buffer);
      
      /* Allocate space                                                 */
      if(cssr == NULL)
      {
         INIT(cssr, CSSR);
         p = cssr;
      }
      else
      {
         ALLOCNEXT(p, CSSR);
      }
      
      /* Check allocation                                               */
      if(p==NULL)
      {
         if(cssr != NULL) FREELIST(cssr,CSSR);
         return(NULL);
      }

      /* V1.1: Initialise all links to zero                             */
      p->link[0] = p->link[1] = p->link[2] = p->link[3] = 
      p->link[4] = p->link[5] = p->link[6] = p->link[7] = 0;
      
      sscanf(buffer,"%d%s%lf%lf%lf%d%d%d%d%d%d%d%d",&p->atnum,
                                                     p->atnam,
                                                    &p->x,
                                                    &p->y,
                                                    &p->z,
                                                    &p->link[0],
                                                    &p->link[1],
                                                    &p->link[2],
                                                    &p->link[3],
                                                    &p->link[4],
                                                    &p->link[5],
                                                    &p->link[6],
                                                    &p->link[7]);
      if(!nocharges)
         sscanf(buffer+73,"%lf",&p->charge);
      else
         p->charge = 0.0;
         
      sscanf(buffer+83,"%d",&p->group);
      
   }
   
   /* Normalise if necessary                                            */
   if(!orthonormal) NormaliseCSSR(cssr,cell,alpha,beta,gamma);
   
   return(cssr);
}

/************************************************************************/
/*>PDB *ReadCSSRasPDB(FILE *fp, int *natom)
   ----------------------------------------
   Input:   FILE     *fp      A pointer to type FILE in which the
                              CSSR file is stored.
   Output:  int      *natom   Number of atoms read.
   Returns: PDB      *pdb     A pointer to the first allocated item of
                              the PDB linked list
   Read a CSSR file into a PDB linked list

   06.09.91 Original
   01.06.92 Documented
   10.06.93 Returns TRUE or FALSE to indicate success
   09.07.93 Changed allocation scheme. Now returns pointer to start of
            list. No need to call init_pdb
   13.07.93 Returns NULL if allocation failed
   27.07.93 Changed I/O to double precision
*/
PDB *ReadCSSRasPDB(FILE  *fp,
                   int   *natom)
{
   char  buffer[160],
         code[8];
   REAL  cell[3],
         alpha, beta, gamma;
   int   orthonormal,
         nocharges,
         i,
         link[8];
   char  title[80],
         name[80];
   PDB   *pdb = NULL,
         *p;
   
   /* Initilialise cell parameters, etc.                                */
   cell[0]  = cell[1]  = cell[2]  = 1.0;
   alpha    = beta     = gamma    = 90.0;
   name[0]  = title[0] = '\0';
   
   *natom   = 0;
   
   /*** Record 1                                                      ***/
   if(!fgets(buffer,160,fp)) return(NULL);
   TERMINATE(buffer);
   
   if(strlen(buffer) >= 28)
   {
      /* We can get a code name                                         */
      strncpy(code,buffer+22,6);
      code[6] = '\0';
   }
   
   if(strlen(buffer) >= 62)
   {
      /* We can get unit cell parameters                                */
      sscanf(buffer+38,"%lf%lf%lf",&cell[0],&cell[1],&cell[2]);
   }
   
   /*** Record 2                                                      ***/
   if(!fgets(buffer,160,fp)) return(NULL);
   TERMINATE(buffer);
   
   if(strlen(buffer) >= 45)
   {
      /* We can get cell angles                                         */
      sscanf(buffer+21,"%lf%lf%lf",&alpha,&beta,&gamma);
   }
   
   /*** Record 3                                                      ***/
   if(!fgets(buffer,160,fp)) return(NULL);
   TERMINATE(buffer);
   
   sscanf(buffer,"%d%d",natom,&orthonormal);
   strcpy(name,buffer+9);
   
   /*** Record 4                                                      ***/
   if(!fgets(buffer,160,fp)) return(NULL);
   TERMINATE(buffer);
   
   sscanf(buffer+4,"%d",&nocharges);
   if(nocharges != 1) nocharges = 0;
   strcpy(title,buffer+7);
   
   /*** Remaining records represent the atoms...                      ***/
   for(i=0; i < *natom; i++)
   {
      if(!fgets(buffer,160,fp)) return(pdb);
      TERMINATE(buffer);
      
      /* Allocate space                                                 */
      if(pdb == NULL)
      {
         INIT(pdb, PDB);
         p = pdb;
      }
      else
      {
         ALLOCNEXT(p, PDB);
      }
      
      /* Check allocation                                               */
      if(p==NULL)
      {
         if(pdb==NULL) FREELIST(pdb, PDB);
         return(NULL);
      }

      sscanf(buffer,"%d%s%lf%lf%lf%d%d%d%d%d%d%d%d",&p->atnum,
                                                     p->atnam,
                                                    &p->x,
                                                    &p->y,
                                                    &p->z,
                                                    &link[0],
                                                    &link[1],
                                                    &link[2],
                                                    &link[3],
                                                    &link[4],
                                                    &link[5],
                                                    &link[6],
                                                    &link[7]);
      strcpy(p->record_type,"ATOM   ");
      padterm(p->atnam,4);
      strcpy(p->resnam,"ATM ");
      strcpy(p->insert," ");
      strcpy(p->chain," ");
      p->occ = 1.0;
      p->resnum = 1;
      p->next = NULL;
      
      if(!nocharges)
         sscanf(buffer+73,"%lf",&p->bval);
      else
         p->bval = 0.0;
   }
   
   /* Normalise if necessary                                            */
   if(!orthonormal) NormalisePDB(pdb,cell,alpha,beta,gamma);
   
   return(pdb);
}

/************************************************************************/
/*>void NormaliseCSSR(CSSR *cssr, REAL cell[3], REAL alpha, 
                      REAL beta, REAL gamma)
   --------------------------------------------------------
   I/O:     CSSR     *cssr    Pointer to CSSR linked list
   Input:   REAL     cell[3]  Unit cell dimensions
                     alpha    Unit cell angles
                     beta
                     gamma
   Convert a CSSR linked list in fractional coordinates to orthonormal

   06.09.91 Original
   01.06.92 Documented
   10.06.93 void return
*/
void NormaliseCSSR(CSSR *cssr,
                   REAL cell[3],
                   REAL alpha,
                   REAL beta,
                   REAL gamma)
{
   int   ncode = 1,           /* Assume a* along X and c* || Z          */
         isw   = 0;           /* Fractional-->Orthonormal               */
   REAL  matrix[3][3],        /* Normalisation matrix                   */
         tempx,tempy,tempz;   /* Used during matrix multiplication      */
   CSSR  *p;
   
   ortho(cell,alpha,beta,gamma,matrix,isw,ncode);
   
   /* Now multiply the coordinates by the matrix                        */
   for(p=cssr;p;NEXT(p))
   {
/**** Original code has this which uses only the RU part of the matrix...
      tempz = p->x * matrix[2][0] +
              p->y * matrix[2][1] +
              p->z * matrix[2][2];
      tempy = p->x * matrix[1][0] +
              p->y * matrix[1][1];
      tempx = p->x * matrix[0][0];
****/
      tempx = p->x * matrix[0][0] +
              p->y * matrix[0][1] +
              p->z * matrix[0][2];
      tempy = p->x * matrix[1][0] +
              p->y * matrix[1][1] +
              p->z * matrix[1][2];
      tempz = p->x * matrix[2][0] +
              p->y * matrix[2][1] +
              p->z * matrix[2][2];

      
      p->x = tempx;
      p->y = tempy;
      p->z = tempz;
   }
}

/************************************************************************/
/*>void NormalisePDB(PDB *pdb, REAL cell[3], REAL alpha,
                     REAL beta, REAL gamma)
   -----------------------------------------------------
   I/O:     PDB      *pdb     Pointer to PDB linked list
   Input:   REAL     cell[3]  Unit cell dimensions
                     alpha    Unit cell angles
                     beta
                     gamma
   Convert a PDB linked list in fractional cooridinates to orthonormal

   06.09.91 Original
   01.06.92 Documented
   10.06.93 void return
*/
void NormalisePDB(PDB   *pdb,
                  REAL  cell[3],
                  REAL  alpha,
                  REAL  beta,
                  REAL  gamma)
{
   int   ncode = 1,           /* Assume a* along X and c* || Z          */
         isw   = 0;           /* Fractional-->Orthonormal               */
   REAL  matrix[3][3],        /* Normalisation matrix                   */
         tempx,tempy,tempz;   /* Used during matrix multiplication      */
   PDB   *p;
   
   ortho(cell,alpha,beta,gamma,matrix,isw,ncode);
   
   /* Now multiply the coordinates by the matrix                        */
   for(p=pdb;p;NEXT(p))
   {
      tempx = p->x * matrix[0][0] +
              p->y * matrix[0][1] +
              p->z * matrix[0][2];
      tempy = p->x * matrix[1][0] +
              p->y * matrix[1][1] +
              p->z * matrix[1][2];
      tempz = p->x * matrix[2][0] +
              p->y * matrix[2][1] +
              p->z * matrix[2][2];
   
      
      p->x = tempx;
      p->y = tempy;
      p->z = tempz;
   }
}

/************************************************************************/
/*>void ortho(REAL cell[3], REAL alpha, REAL beta, REAL gamma,
              REAL amatrx[3][3], int isw, int ncode)
   -----------------------------------------------------------
   Input:   REAL     cell[3]  Unit cell dimensions
                     alpha    Unit cell angles
                     beta
                     gamma
   Output:  REAL     amatrx[3][3]   Returned conversion matrix
   Input:   int      isw      0: Frac-->Ortho,  1: Ortho-->Frac
            int      ncode    Orientation of reciprocal axes wrt true axes.
                              1: a || xo,          c* || zo
                              2: b || xo,          a* || zo
                              3: c || xo,          b* || zo
                              4: hex a & b || xo,  c* || zo
                              5: a* || xo,         c  || zo

   Function to calculate a matrix which will convert between fractional and
   orthonormal coordinates given unit cell dimensions and angles.
   ncode defines the orientation of the reciprocal axes wrt the real axes.
   
   See Rollett `Computing Methods in Crystallography' p.23

   06.09.91 Original
   01.06.92 Documented
   10.06.93 void return
*/
void ortho(REAL  cell[3],        /* Cell dimensions                     */
           REAL  alpha,          /* Cell angles                         */
           REAL  beta,
           REAL  gamma,
           REAL  amatrx[3][3],   /* Returned conversion matrix          */
           int   isw,            /* 0: Frac-->Ortho,  1: Ortho-->Frac   */
           int   ncode)          /* Orientation of reciprocal axes      */
{
   REAL  ro[3][3],
         rf[3][3],
         conv  = PI/180.0,
         a     = cell[0],
         b     = cell[1],
         c     = cell[2],
         sina,    cosa,
         sinb,    cosb,
         sing,    cosg,
         sinas,   cosas,
         sinbs,   cosbs,
         sings,   cosgs;
   int   i, j;
         
   /* Convert angles to rad's                                           */
   alpha *= conv;
   beta  *= conv;
   gamma *= conv;
   
   /* Check code                                                        */
   if(ncode==0) ncode = 1;
        
        
   /* Get sin's and cos's                                               */
   sina  = (REAL)sin((double)alpha);
   cosa  = (REAL)cos((double)alpha);
   sinb  = (REAL)sin((double)beta);
   cosb  = (REAL)cos((double)beta);
   sing  = (REAL)sin((double)gamma);
   cosg  = (REAL)cos((double)gamma);
   cosbs = (cosa*cosg-cosb)/(sina*sing);
   sinbs = sqrt(1.0-cosbs*cosbs);
   cosas = (cosg*cosb-cosa)/(sinb*sing);
   sinas = sqrt(1.0-cosas*cosas);
   cosgs = (cosa*cosb-cosg)/(sina*sinb);
   sings = sqrt(1.0-cosgs*cosgs);
   
   for(i=0;i<3;i++)
      for(j=0;j<3;j++)
         ro[i][j] = 0.0;
         
   switch(ncode)
   {
   case 1:
      /* Brookhaven convention                                          */
      ro[0][0] = a;
      ro[0][1] = b*cosg;
      ro[0][2] = c*cosb;
      ro[1][0] = 0.0;
      ro[1][1] = b*sing;
      ro[1][2] = -c*sinb*cosas;
      ro[2][0] = 0.0;
      ro[2][1] = 0.0;
      ro[2][2] = c*sinb*sinas;
      break;
   case 2:
      ro[0][0] = a*cosg;
      ro[0][1] = b;
      ro[0][2] = c*cosa;
      ro[1][0] = -a*sing*cosbs;
      ro[1][1] = 0.0;
      ro[1][2] = c*sina;
      ro[2][0] = a*sing*sinbs;
      ro[2][1] = 0.0;
      ro[2][2] = 0.0;
      break;
   case 3:
      ro[1][0] = a*sinb;
      ro[1][1] = -b*sina*cosgs;
      ro[1][2] = 0.0;
      ro[2][0] = 0.0;
      ro[2][1] = b*sina*sings;
      ro[2][2] = 0.0;
      ro[0][0] = a*cosb;
      ro[0][1] = b*cosa;
      ro[0][2] = c;
      break;
   case 4:
      /* Gerard Bricogne's orthogonal System for P3:
         xo =  a + b
         yo = -a + b
      */
      ro[0][0] = a/2.0;
      ro[0][1] = a/2.0;
      ro[0][2] = 0.0;
      ro[1][0] = -a*sing;
      ro[1][1] = a*sing;
      ro[1][2] = 0.0;
      ro[2][0] = 0.0;
      ro[2][1] = 0.0;
      ro[2][2] = c;
      break;
   case 5:
      /* Rollet matrix as used in Cambridge                             */
      ro[0][0] = a*sinb*sings;
      ro[0][1] = 0.0;
      ro[0][2] = 0.0;
      ro[1][0] = -a*sinb*cosgs;
      ro[1][1] = b*sina;
      ro[1][2] = 0.0;
      ro[2][0] = a*cosb;
      ro[2][1] = b*cosa;
      ro[2][2] = c;
      break;
   }

   invert33(ro,rf);
   for(i=0;i<3;i++)
   {
      for(j=0;j<3;j++)
      {
         if(isw)
            amatrx[j][i] = rf[j][i];
         else
            amatrx[j][i] = ro[j][i];
      }
   }

}

