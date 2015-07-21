/*************************************************************************

   Program:    
   File:       HAddPDB.c
   
   Version:    V2.16R
   Date:       24.01.06
   Function:   Add hydrogens to a PDB linked list
   
   Copyright:  (c) SciTech Software 1990-2006
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
   Routine to add hydrogens to a protein linked list of type PDB.
   The routine allocates space for the new atoms and inserts them
   into the list at the appropriate positions within the residues.

   N.B. Because of the variation in styles required for N-terminal Hs,
   this routine adds *no* hydrogens to the N-terminal residue. FixPDB()
   may be used for this purpose if CHARMM/CONGEN style hydrogens are
   required.

**************************************************************************

   Usage:
   ======
   nhyd = HAddPDB(fp,pdb)
   Input:         FILE   *fp        File containing proton generation
                                    parameters.
   Input/Output:  PDB    *pdb       Linked list of protein structure.
   Returns:       int               Number of hydrogens added.

   The globally defined structure gHaddInfo gives information on the 
   number of each hydrogen type created. This structure is defined as
   follows:
   typdef struct
   {
      int   Total,      Total hydrogens
            T1,         Type 1 C-H's
            T2,         Type 2 C-H2's
            T3,         Type 3 C-H3's
            T4,         Type 4 sp2 C-H's,>N-H
            T5;         Type 5 O-H's =N-H's
   }  HADDINFO;
   To examine these values in your code, reference the structure as:
   extern HADDINFO gHaddInfo;

**************************************************************************

   Revision History:
   =================
   V2.0  16.05.90 AddH is changed to insert each set of atoms for each 
                  PGP, on the fly, rather than building a complete list 
                  of hydrogens and then merging the two lists. This allows
                  us to get round the problem of missing atoms, since 
                  there will be no merging error.

   V2.1  24.05.90 Returns the number fo hydrogens added. Also fixes bug 
                  relating to number of type 2 and type 3 H's added.
                  Doesn't work under UNIX!

   V2.2  15.07.91 A few bits of tidying up:
                  >  Now uses macros.h rather than defining macros itself. 
                  >  Arrays now changed so should fix alignment problems 
                     under UNIX. 
                  >  Now correctly checks return from forscanf() and no 
                     longer reads characters to check EOF itself.
                  >  Improves treatment of NTER residues where it now 
                     generates H's on the first true residue. Residues 
                     labelled NTER will be ignored. Calling FixNterH() 
                     will move the H coords into the NTER residue if 
                     required. 
                  >  Fixes reading of PGP files with blank lines 
                  >  Bug fix for type 3's
                  Currently untested under UNIX.

   V2.3  27.07.93 Changed to use fsscanf() and I/O precision is double
   V2.4  08.03.94 Changed static variable names and added casts on
                  maths functions. Added Dummy atom handling in makeh().
                  A few other bits of tidying.
   V2.5  23.08.94 Added OpenPGPFile() routine
   V2.6  01.09.94 Fixed bug in some compilers in ReadPGP()
   V2.7  26.01.96 Wasn't handling insert codes in PDB.
   V2.8  24.05.99 Fixed two memory leaks in makeh()
                  Also skips HETATMs in GenH()
   V2.9  30.05.02 Changed PDB field from 'junk' to 'record_type'
   V2.10 05.12.02 Correctly sets the atnam_raw field
   V2.11 27.03.03 Fixed severe memory leak in AddH()
   V2.12 03.06.05 Added altpos
   V2.13 28.07.05 Added conditionals for msdos and Mac OS/X
   V2.14 28.11.05 No longer exits if previous C is missing
   V2.15 24.01.06 Fixed error message in GenH() which could try to print 
                  from NULL pointer

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
#include "fsscanf.h"
#include "macros.h"
#include "general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXTYPE 200        /* Max number of H definitions in PGP file   */
#define DATAENV "DATADIR"         /* Unix environment variable for data */
#define DATADIR "AMDATA:"         /* VMS/AMigaDOS assign for data       */
#define EXPLPGP "Explicit.pgp"    /* The PGP filename                   */
#define ALLHPGP "AllH.pgp"        /* The PGP filename                   */
/* #define SCREEN_INFO */  /* Error messages output via screen()        */

/************************************************************************/
/* Define-dependent includes
*/
#ifdef SCREEN_INFO
#   include "WindIO.h"
#endif

/************************************************************************/
/* Globals local to this file
*/
static char    sGRes[MAXTYPE][8],
               sGAtom[MAXTYPE][8][8],
               sGRName[8],
               sGHName[8][8],
               sGNat[16][8],
               sIns;
static int     sHType[MAXTYPE],
               sNpgp,
               sNo,
               sKMax,
               sNType1, sNType2, sNType3, sNType4, sNType5;
static REAL    sGR[MAXTYPE],
               sAlpha[MAXTYPE],
               sBeta[MAXTYPE],
               sGX[16], sGY[16], sGZ[16],
               sFac;

/************************************************************************/
/* Globals which are externally visible                                  
*/
HADDINFO gHaddInfo;

/************************************************************************/
/* Prototypes for static function
*/
static int GenH(PDB *pdb, BOOL *err_flag);
static PDB *makeh(int HType, REAL BondLen, REAL alpha, REAL beta, 
                  BOOL firstres);
static BOOL AddH(PDB *hlist, PDB **position, int HType);
static void SetRawAtnam(char *out, char *in);
static PDB  *StripDummyH(PDB *pdb, int *nhyd);

/************************************************************************/
/*>int HAddPDB(FILE *fp, PDB  *pdb)
   --------------------------------
   Input:   FILE     *fp       File pointer to PGP file
   I/O:     PDB      *pdb      PDB Linked list to which Hs are added
   Returns: int                Number of Hs added. 0 if error.
   Globals: HADDINFO gHaddInfo Information on Hs added

   This routine adds hydrogens to a PDB linked list. Performs all
   necessary functions.

   N.B. Because of the variation in styles required for N-terminal Hs,
   this routine adds *no* hydrogens to the N-terminal residue. FixPDB()
   may be used for this purpose if CHARMM/CONGEN style hydrogens are
   required.

   16.05.90 Original    By: ACRM
   04.01.94 Changed check on return=NULL to 0
   08.03.94 Only reads PGP on first call. err_flag changed to BOOL.
   28.11.05 Removes any dummy hydrogens added because there were
            missing atoms
*/
int HAddPDB(FILE *fp, PDB  *pdb)
{
   PDB            *p;
   int            nhydrogens,
                  atomcount = 1;
   BOOL           err_flag  = FALSE;
   static BOOL    FirstCall = TRUE;
   

   if(FirstCall)
   {
      /* Read the parameter file                                        */
      FirstCall = FALSE;
      ReadPGP(fp);
   }

   /* Generate the hydrogens                                            */
   if((nhydrogens=GenH(pdb,&err_flag))==0)
      return(0);

/* ACRM+++ 28.11.05                                                     */
   /* Remove dummy hydrogens (where atoms are missing)                  */
   pdb = StripDummyH(pdb, &nhydrogens);
/* ACRM=== 28.11.05                                                     */
   
   /* Renumber atoms in PDB linked list                                 */
   for(p=pdb;p;NEXT(p)) p->atnum=atomcount++;
   return(nhydrogens);
}

/************************************************************************/
/*>int ReadPGP(FILE *fp)
   ---------------------
   Input:   FILE *fp  Pointer to PGP file.
   Returns: int       Number of parameters read.

   Read a proton generation parameter file. All data are placed in static
   arrays used by the hydrogen adding routines.
   Normally this routine is only called by the first call to HAddPDB().
   It is only necessary to call this routine explicitly if the PGP file
   is changed between calls to HAddPDB() and thus needs re-reading.

   16.05.90 Original    By: ACRM
   27.07.93 Changed to use fsscanf()
   01.03.94 Changed static variable names
   01.09.94 Moved n++ out of the fsscanf()
*/
int ReadPGP(FILE *fp)
{
   char  buffer[160];
   int   n=0;
   
   while(fgets(buffer,159,fp))
   {
      n++;
      
      fsscanf(buffer,
              "%4s%4s%1x%4s%1x%4s%1x%4s%1x%4s%1x%4s%1x%1d%10lf%10lf%10lf",
              sGRes[n],              
              sGAtom[n][1],
              sGAtom[n][2],
              sGAtom[n][3],
              sGAtom[n][4],
              sGAtom[n][5],
              sGAtom[n][6],
              &(sHType[n]),
              &(sGR[n]),
              &(sAlpha[n]),
              &(sBeta[n]));

#ifdef DEBUG_READ
      printf("%4s %4s %4s %4s %4s %4s %4s %1d %8.3f %8.3f %8.3f\n",
             sGRes[n],sGAtom[n][1],sGAtom[n][2],sGAtom[n][3],sGAtom[n][4],
             sGAtom[n][5],sGAtom[n][6],sHType[n],sGR[n],
             sAlpha[n],sBeta[n]);
#endif

      if(sHType[n] != 0)
      {
         sAlpha[n] *= (PI/180.0);    
         sBeta[n]  *= (PI/180.0);
      }
   }  /* End of file                                                    */

   sNpgp = n;
   return(sNpgp);
}

/************************************************************************/
/*>static int GenH(PDB *pdb, BOOL *err_flag)
   -----------------------------------------
   I/O:     PDB  *pdb      PDB Linked to which Hs are added
            BOOL *err_flag Error flag
   Returns: int            Number of hydrogens added (0 if error)

   Does the actual work of generating a set of hydrogens

   16.05.90 Original    By: ACRM
   01.03.94 Changed static variable names
   08.03.94 Changed err_flag to BOOL
   26.01.96 Added check on insert code as well as resnum
   28.11.05 Modified such that a missing preceeding C no longer causes
            the routine to exit
   24.01.06 Fixed error message which could try to print from NULL 
            pointer
*/
static int GenH(PDB *pdb, BOOL *err_flag)
{
   BOOL  firstres;
   char  *bl       = "    ",
         *co       = "CO  ",
         *c        = "C   ";
   int   k, n, m, j, ittot;
   PDB   *p,*q;
   PDB   *position[16],*hlist;
#ifdef SCREEN_INFO
   char  buffer[160];
#endif
    
   for(j=0;j<16;position[j++]=NULL);

   sFac     = 0.5 * (REAL)sqrt((double)3.0);
    
   firstres = TRUE;
    
   /* Main loop                                                         */
   sNType1=0; sNType2=0; sNType3=0; sNType4=0; sNType5=0;

   /* For the first residue, we don't have info for the previous C
      so set sGNat[1] (the atom list for this residue) to a blank
   */
   strcpy(sGNat[1],bl);
   
   position[1]=NULL;

   /* For each atom in the PDB file                                     */
   for(p=pdb;p;)
   {
      if(!strncmp(p->record_type,"HETATM",6))      /* 24.05.99          */
      {
         NEXT(p);
         continue;
      }
      

      /* Don't do anything with NTER residues except reset firstres     */
      if(!strncmp(p->resnam,"NTER",4))
      {
         firstres = TRUE;
         NEXT(p);
         continue;
      }

      k=1;  /* This is one as this position is reserved for the CO of
               the previous residue 
            */

      /* Copy this residue into our global work arrays.
         sGRName         is the residue name
         sGX[],sGY[],sGZ[] are the coordinates of the atoms
         sGNat[]         are the atom names
      */
      do
      {
         k++;
         position[k] = p;

         sNo=p->resnum;
         sIns=p->insert[0];                   /* 26.01.96               */
         strcpy(sGRName,p->resnam);
#ifdef DEBUG
         printf("Group name is: %s\n",sGRName);
#endif
         strcpy(sGNat[k],p->atnam);
#ifdef DEBUG
         printf("Atom %d name: %s\n",k,sGNat[k]);
#endif
         sGX[k] = p->x;
         sGY[k] = p->y;
         sGZ[k] = p->z;
         NEXT(p);
         if(!p) break;
      } while((p != NULL)        && 
              (p->resnum == sNo) && 
              (p->insert[0] == sIns));       /* 26.01.96                */

      /* sKMax is used to store the number of atoms in this residue     */
      sKMax = k;

      /* Go through each of the PGP's until we find this residue type   */
      for(n=1;n<=sNpgp;n++)
      {
#ifdef DEBUG
         printf("sGRes(%d)=%s\n",n,sGRes[n]);
#endif
         if(strncmp(sGRName,sGRes[n],4)) continue;
#ifdef DEBUG
         printf("Entry found for this res. in PGP\n");
#endif

         /* Having found the residue type, copy the associated PGP atom
            list into the sGHName[] work array
         */
         for(m=1;m<=6;m++) strcpy(sGHName[m],sGAtom[n][m]);
         
         /* Now generate the hydrogen(s) associated with this PGP       */
         if((hlist = makeh(sHType[n],sGR[n],sAlpha[n],sBeta[n],firstres))
            !=NULL)
         {
            /* And add it into the list, updating p to point to the new
               end of this residue
            */
            if(!AddH(hlist,position,sHType[n])) return(0);
         }
         if(*err_flag) return(0);
      }
      
      /* If this is the first residue then handle it as NTER            */
      if(firstres)
      {
#ifdef DEBUG
         printf("Rechecking for NTER...\n");
#endif
         for(n=1;n<=sNpgp;n++)
         {
#ifdef DEBUG
            printf("sGRes(%d)=%s\n",n,sGRes[n]);
#endif
            if(strncmp("NTER",sGRes[n],4)) continue;
#ifdef DEBUG
            printf("Entry found for NTER in PGP.\n");
#endif

            /* Having found the residue type, copy the associated PGP atom
               list into the sGHName[] work array                        
            */
            for(m=1;m<=6;m++) strcpy(sGHName[m],sGAtom[n][m]);
         
            /* Now generate the hydrogen associated with this PGP       */
            if((hlist = makeh(sHType[n],sGR[n],sAlpha[n],
                              sBeta[n],firstres))!=NULL)
            {
               /* And add it into the list, updating p to point to the new
                  end of this residue 
               */
               if(!AddH(hlist,position,sHType[n])) return(0);
            }
            if(*err_flag) return(0);
         }
      }

      /* Next amino acid
         Set up pointer for carbonyl of previous residue unless this 
         residue is a CTER when we set firstres to TRUE
      */
      
      if(strncmp(sGRName,"CTER",4))
      {
         firstres = FALSE;
         for(j=sKMax-1;j>0;j--) if(!strncmp(sGNat[j],c,4))break;
         if(j==0)
         {
#ifdef SCREEN_INFO
            sprintf(buffer,"\nError==> genh() found no carbonyl carbon \
in residue %d\n\n",p->resnum);
            screen(buffer);
#endif
/* ACRM--- 25.11.05
            *err_flag=TRUE;
            return(0);
*/
/* ACRM+++ 28.11.05                                                     */
/* ACRM 24.01.06                                                        */
            if(p!=NULL)
            {
               fprintf(stderr,"Warning=> genh() found no carbonyl carbon \
preceeding residue %c%d%c\n", p->chain[0], p->resnum, p->insert[0]);
            }
            else
            {
               fprintf(stderr,"Warning=> genh() found no carbonyl carbon \
preceeding the last residue\n");
            }
            sGX[1]=9999.0;
            sGY[1]=9999.0;
            sGZ[1]=9999.0;
            
/* ACRM=== 28.11.05                                                     */
         }
         else
         {
            sGX[1]=sGX[j];
            sGY[1]=sGY[j];
            sGZ[1]=sGZ[j];
         }
         strcpy(sGNat[1],co);
         q=position[j];
            
         for(k=0;k<16;position[k++]=NULL);
         position[1]=q;
      }
      else
      {
         firstres = TRUE;
      }
      
   }  /* Go back to the next atom/residue                               */
   
   ittot=sNType1+sNType2+sNType3+sNType4+sNType5;
   
   gHaddInfo.Total = ittot;
   gHaddInfo.T1    = sNType1;
   gHaddInfo.T2    = sNType2;
   gHaddInfo.T3    = sNType3;
   gHaddInfo.T4    = sNType4;
   gHaddInfo.T5    = sNType5;

   return(ittot);
}

/************************************************************************/
/*>static PDB *makeh(int HType, REAL BondLen, REAL alpha, REAL beta, 
                     BOOL firstres)
   -----------------------------------------------------------------
   Input:   int  HType    Hydrogen type number
            REAL BondLen  Length of the bond to the added hydrogen
            REAL alpha    The first angle defining the H coordinate
            REAL beta     The second angle defining the H coordinate
            REAL firstres If set, don't add a planar hydrogen as it's 
                          the first residue of a chain
   Returns: PDB *         Pointer to linked list of created hydrogens.
                          NULL if memory allocation fails or this is the 
                          Nter N where we don't require a planar H or
                          all atoms have been done.

   Generate a set of hydrogen coordinates. The antecedent atoms are set
   in static external variables. A linked list containing the 
   appropriate number of Hs is returned.

   16.05.90 Original    By: ACRM
   08.03.94 Added code to handle dummy atom positions (All occurences
            of variable `Dummy'). Made assigned character strings static.
   26.01.96 Now stores insert codes into hlist
   24.05.99 Fixed two memory leaks
   05.12.02 Added setting of atnam_raw
   03.06.05 Added setting of altpos
*/
static PDB *makeh(int HType, REAL BondLen, REAL alpha, REAL beta, 
                  BOOL firstres)
{ 
   static char    *nt = "NT  ",
                  *n  = "N   ";
   REAL           x1,y1,z1,x2,y2,z2,x3,y3,z3,
                  x4,y4,z4,x5,y5,z5,x6,y6,z6,
                  x21,y21,z21,r21,
                  x21p,y21p,z21p,r21p,
                  x23,y23,z23,r23,
                  xp23,yp23,zp23,rp23,
                  x24,y24,z24,r24,
                  x32,y32,z32,
                  xv25,yv25,zv25,rv25,
                  cosa,sina,
                  cosb,sinb,
                  cosax,cosay,cosaz,
                  xa,ya,za,xb,yb,zb,
                  xab,yab,zab,rab,
                  xmin,ymin,zmin,
                  xapb,yapb,zapb,rapb,
                  xplus,yplus,zplus,
                  xp,yp,zp,
                  xs,ys,zs,
                  xh,yh,zh,
                  xv,yv,zv,
                  scalpr;
   int            kount = 0,
                  num_ant,k,jj;
   unsigned short nt_point;
   BOOL           ok,
                  Dummy;
   PDB            *hlist,
                  *hlist_start = NULL;

   Dummy = FALSE;

   INIT(hlist_start, PDB);
   if(hlist_start == NULL) return(NULL);
   
   CLEAR_PDB(hlist_start);
   hlist = hlist_start;

   nt_point = 0;
   if(HType==1) num_ant=4; else num_ant=3;

   /* Don't add a planar H to the Nter N.                               */
   if(firstres && HType==4 && !strncmp(sGHName[2],n,4)) 
   {
      /* 24.05.99 Fixed memory leak                                     */
      free(hlist_start);
      return(NULL);
   }

   /* Work through the atoms in this residue (sGNat[]) and compare them 
      with the first 4 atoms in the PGP atom list in sGHName[], storing 
      the associated coordinates
   */
   for(k=1; k<=sKMax; k++)
   {
      if(nt_point) 
      {
         strcpy(sGNat[nt_point],"N   ");
         nt_point=0;
      }

      if(!strncmp(sGHName[1],sGNat[k],4))
      {
         if(!strncmp(sGNat[k],"NT  ",4)) nt_point=k;
         kount++;
         x1=sGX[k];
         y1=sGY[k];
         z1=sGZ[k];

         /* Check for dummy atom                                        */
         if(x1 > (REAL)9998.0 && 
            y1 > (REAL)9998.0 && 
            z1 > (REAL)9998.0)
            Dummy = TRUE;
      }
      else if(!strncmp(sGHName[2],sGNat[k],4))
      {
         if(!strncmp(sGNat[k],"NT  ",4)) nt_point=k;
         kount++;
         x2=sGX[k];
         y2=sGY[k];
         z2=sGZ[k];

         /* Check for dummy atom                                        */
         if(x2 > (REAL)9998.0 && 
            y2 > (REAL)9998.0 && 
            z2 > (REAL)9998.0)
            Dummy = TRUE;
      }
      else if(!strncmp(sGHName[3],sGNat[k],4))
      {
         if(!strncmp(sGNat[k],"NT  ",4)) nt_point=k;
         kount++;
         x3=sGX[k];
         y3=sGY[k];
         z3=sGZ[k];

         /* Check for dummy atom                                        */
         if(x3 > (REAL)9998.0 && 
            y3 > (REAL)9998.0 && 
            z3 > (REAL)9998.0)
            Dummy = TRUE;
      }
      else if(!strncmp(sGHName[4],sGNat[k],4))
      {
         if(!strncmp(sGNat[k],"NT  ",4)) nt_point=k;
         kount++;
         x4=sGX[k];
         y4=sGY[k];
         z4=sGZ[k];

         /* Check for dummy atom                                        */
         if(x4 > (REAL)9998.0 && 
            y4 > (REAL)9998.0 && 
            z4 > (REAL)9998.0)
            Dummy = TRUE;
      }
   }  /* End of k loop around this residue                              */

   /* Check we found all the atoms we need for this PGP                 */
   if(kount != num_ant)
   {
      ok = FALSE;
      if(firstres)
      {
         /* Check it's not the missing N in the first residue           */
         for(jj=1;jj<=4;jj++) if(!strncmp(sGHName[jj],n,4)) ok = TRUE;
      }
      else
      {
         /* Check it's not just the NT                                  */
         for(jj=1;jj<=4;jj++) if(!strncmp(sGHName[jj],nt,4)) ok = TRUE;
      }
#ifdef SCREEN_INFO
      if(!ok)
      {
         char buffer[160];
         
         sprintf(buffer,"Error==> makeh() unable to find all atoms \
required by PGP parameter for %3s %5d%c\n",sGRName,sNo,sIns);
         screen(buffer);
         screen("Atoms required by PGP\n");
         screen("SGHNAME: ");
         for(jj=1;jj<=4;jj++)
         {
            sprintf(buffer," %4s",sGHName[jj]);
            screen(buffer);
         }
         screen("\n");
         screen("Atoms in current residue\n");
         screen("SGNAT  : ");
         for(jj=1;jj<=sKMax;jj++)
         {
            sprintf(buffer," %4s",sGNat[jj]);
            screen(buffer);
         }
         screen("\n");
      }
#endif
      /* 24.05.99 Fixed memory leak                                     */
      FREELIST(hlist_start, PDB);
      return(NULL);
   }
    
   x21=x2-x1;
   y21=y2-y1;
   z21=z2-z1;
   r21=(REAL)sqrt((double)(x21*x21 + y21*y21 + z21*z21));

   x23=x2-x3;
   y23=y2-y3;
   z23=z2-z3;
   r23=(REAL)sqrt((double)(x23*x23 + y23*y23 + z23*z23));
    
   if(HType == 1)
   {
      /* HTYPE 1: Generation of 1 tetrahedral H
         --------------------------------------
      */
    
      if(Dummy)
      {
         x5 = y5 = z5 = (REAL)9999.0;
      }
      else
      {
         x24=x2-x4;
         y24=y2-y4;
         z24=z2-z4;
         r24=(REAL)sqrt((double)(x24*x24 + y24*y24 + z24*z24));
         xv25=x21/r21+x24/r24+x23/r23;
         yv25=y21/r21+y24/r24+y23/r23;
         zv25=z21/r21+z24/r24+z23/r23;
         rv25=(REAL)sqrt((double)(xv25*xv25 + yv25*yv25 + zv25*zv25));
         x5=x2+BondLen*xv25/rv25;
         y5=y2+BondLen*yv25/rv25;
         z5=z2+BondLen*zv25/rv25;
      }

      hlist->resnum=sNo;
      hlist->insert[0]=sIns;              /* 26.01.96                   */
      padterm(hlist->insert,4);
      strcpy(hlist->atnam,sGHName[5]);
      SetRawAtnam(hlist->atnam_raw, sGHName[5]);   /* 05.12.02          */
      hlist->altpos = ' ';                         /* 03.06.05          */
      hlist->x=x5;
      hlist->y=y5;
      hlist->z=z5;
      ALLOCNEXT(hlist,PDB);
      if(hlist == NULL)
      {
         FREELIST(hlist_start, PDB);
         return(NULL);
      }
      CLEAR_PDB(hlist);

#ifdef DEBUG
      printf("makeh() Type 1 allocated hlist = %d\n",(int)hlist);
#endif

      sNType1++;
   }         /* End of HTYPE 1                                         */
   else      /* All types other than HTYPE 1                           */
   {
      cosa=(REAL)cos((double)alpha);
      sina=(REAL)sin((double)alpha);
      switch(HType)
      {
      case 2:  
         /* HTYPE 2: Generation of 2 tetrahedral H's
            ----------------------------------------
         */
         if(Dummy)
         {
            x4 = y4 = z4 = (REAL)9999.0;
            x5 = y5 = z5 = (REAL)9999.0;
         }
         else
         {
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
            x4=x2+BondLen*(cosa*xplus+sina*xs);
            y4=y2+BondLen*(cosa*yplus+sina*ys);
            z4=z2+BondLen*(cosa*zplus+sina*zs);
            x5=x2+BondLen*(cosa*xplus-sina*xs);
            y5=y2+BondLen*(cosa*yplus-sina*ys);
            z5=z2+BondLen*(cosa*zplus-sina*zs);
         }

         hlist->resnum=sNo;
         hlist->insert[0]=sIns;           /* 26.01.96                   */
         padterm(hlist->insert,4);
         strcpy(hlist->atnam,sGHName[4]);
         SetRawAtnam(hlist->atnam_raw, sGHName[4]);   /* 05.12.02       */
         hlist->altpos = ' ';                         /* 03.06.05       */
         hlist->x=x4;
         hlist->y=y4;
         hlist->z=z4;
         ALLOCNEXT(hlist,PDB);
         if(hlist == NULL)
         {
            FREELIST(hlist_start, PDB);
            return(NULL);
         }
         CLEAR_PDB(hlist);

#ifdef DEBUG
         printf("makeh() Type 2a allocated hlist = %d\n",(int)hlist);
#endif

         hlist->resnum=sNo;
         hlist->insert[0]=sIns;           /* 26.01.96                   */
         padterm(hlist->insert,4);
         strcpy(hlist->atnam,sGHName[5]);
         SetRawAtnam(hlist->atnam_raw, sGHName[5]);   /* 05.12.02       */
         hlist->altpos = ' ';                         /* 03.06.05       */
         hlist->x=x5;
         hlist->y=y5;
         hlist->z=z5;
         ALLOCNEXT(hlist,PDB);
         if(hlist == NULL)
         {
            FREELIST(hlist_start, PDB);
            return(NULL);
         }
         CLEAR_PDB(hlist);

#ifdef DEBUG
         printf("makeh() Type 2b allocated hlist = %d\n",(int)hlist);
#endif

         sNType2+=2;
         break;

/* Initialisation for both these cases is the same                      */
      case 3:
      case 5:
         if(!Dummy)
         {
            /* Don't bother with all this lot if we have a dummy atom   */
            x32=x3-x2;
            y32=y3-y2;
            z32=z3-z2;
            xh=x32/r23;
            yh=y32/r23;
            zh=z32/r23;
            scalpr=(x21*x32+y21*y32+z21*z32)/r23;
            xp=scalpr*xh;
            yp=scalpr*yh;
            zp=scalpr*zh;
            x21p=x21-xp;
            y21p=y21-yp;
            z21p=z21-zp;
            r21p=(REAL)sqrt((double)(x21p*x21p+y21p*y21p+z21p*z21p));
            xv=x21p/r21p;
            yv=y21p/r21p;
            zv=z21p/r21p;
            xs=yh*zv-zh*yv;
            ys=zh*xv-xh*zv;
            zs=xh*yv-yh*xv;
            cosax=cosa*xh;
            cosay=cosa*yh;
            cosaz=cosa*zh;
         }

         if(HType==3)
         {
            /* HTYPE 3: Generation of 3 tetrahedral H's
               ----------------------------------------
            */
            if(Dummy)
            {
               x4 = y4 = z4 = (REAL)9999.0;
               x5 = y5 = z5 = (REAL)9999.0;
               x6 = y6 = z6 = (REAL)9999.0;
            }
            else
            {
               x4=x3+BondLen*(cosax+sina*xv);
               y4=y3+BondLen*(cosay+sina*yv);
               z4=z3+BondLen*(cosaz+sina*zv);
               
               /* V2.2: Bug fix here: xy, ys, zs; not xs all the time!  */
               x5=x3+BondLen*(cosax+sina*(sFac*xs-0.5*xv));
               y5=y3+BondLen*(cosay+sina*(sFac*ys-0.5*yv));
               z5=z3+BondLen*(cosaz+sina*(sFac*zs-0.5*zv));
               x6=x3+BondLen*(cosax+sina*(-sFac*xs-0.5*xv));
               y6=y3+BondLen*(cosay+sina*(-sFac*ys-0.5*yv));
               z6=z3+BondLen*(cosaz+sina*(-sFac*zs-0.5*zv));
            }

            hlist->resnum=sNo;
            hlist->insert[0]=sIns;        /* 26.01.96                   */
            padterm(hlist->insert,4);
            strcpy(hlist->atnam,sGHName[4]);
            SetRawAtnam(hlist->atnam_raw, sGHName[4]);   /* 05.12.02    */
            hlist->altpos = ' ';                         /* 03.06.05    */
            hlist->x=x4;
            hlist->y=y4;
            hlist->z=z4;
            ALLOCNEXT(hlist,PDB);
            if(hlist == NULL)
            {
               FREELIST(hlist_start, PDB);
               return(NULL);
            }
            CLEAR_PDB(hlist);

#ifdef DEBUG
            printf("makeh() Type 3a allocated hlist = %d\n",(int)hlist);
#endif

            hlist->resnum=sNo;
            hlist->insert[0]=sIns;        /* 26.01.96                   */
            padterm(hlist->insert,4);
            strcpy(hlist->atnam,sGHName[5]);
            SetRawAtnam(hlist->atnam_raw, sGHName[5]);   /* 05.12.02    */
            hlist->altpos = ' ';                         /* 03.06.05    */
            hlist->x=x5;
            hlist->y=y5;
            hlist->z=z5;
            ALLOCNEXT(hlist,PDB);
            if(hlist == NULL)
            {
               FREELIST(hlist_start, PDB);
               return(NULL);
            }
            CLEAR_PDB(hlist);

#ifdef DEBUG
            printf("makeh() Type 3b allocated hlist = %d\n",(int)hlist);
#endif

            hlist->resnum=sNo;
            hlist->insert[0]=sIns;        /* 26.01.96                   */
            padterm(hlist->insert,4);
            strcpy(hlist->atnam,sGHName[6]);
            SetRawAtnam(hlist->atnam_raw, sGHName[6]);   /* 05.12.02    */
            hlist->altpos = ' ';                         /* 03.06.05    */
            hlist->x=x6;
            hlist->y=y6;
            hlist->z=z6;
            ALLOCNEXT(hlist,PDB);
            if(hlist == NULL)
            {
               FREELIST(hlist_start, PDB);
               return(NULL);
            }
            CLEAR_PDB(hlist);

#ifdef DEBUG
            printf("makeh() Type 3c allocated hlist = %d\n",(int)hlist);
#endif

            sNType3+=3;
         }
         else if(HType==5)
         {
            /* HTYPE 5: Generation of 1 H where an angle is specified
               ------------------------------------------------------
            */
            if(Dummy)
            {
               x4 = y4 = z4 = (REAL)9999.0;
            }
            else
            {
               cosb=(REAL)cos((double)beta);
               sinb=(REAL)sin((double)beta);
               x4=x3+BondLen*(cosax+sina*(cosb*xv+sinb*xs));
               y4=y3+BondLen*(cosay+sina*(cosb*yv+sinb*ys));
               z4=z3+BondLen*(cosaz+sina*(cosb*zv+sinb*zs));
            }

            hlist->resnum=sNo; 
            hlist->insert[0]=sIns;        /* 26.01.96                   */
            padterm(hlist->insert,4);
            strcpy(hlist->atnam,sGHName[4]);
            SetRawAtnam(hlist->atnam_raw, sGHName[4]);   /* 05.12.02    */
            hlist->altpos = ' ';                         /* 03.06.05    */
            hlist->x=x4;
            hlist->y=y4;
            hlist->z=z4;
            ALLOCNEXT(hlist,PDB);
            if(hlist == NULL)
            {
               FREELIST(hlist_start, PDB);
               return(NULL);
            }
            CLEAR_PDB(hlist);

#ifdef DEBUG
            printf("makeh() Type 5 allocated hlist = %d\n",(int)hlist);
#endif

            sNType5++;
         }
         break;

      case 4:  
         /* HTYPE 4: Generation of 1 sp2 H
            ------------------------------
         */
         if(Dummy)
         {
            x4 = y4 = z4 = (REAL)9999.0;
         }
         else
         {
            x32=x3-x2;
            y32=y3-y2;
            z32=z3-z2;
            scalpr=(x21*x32+y21*y32+z21*z32)/r21;
            xh=x21/r21;
            yh=y21/r21;
            zh=z21/r21;
            xp=scalpr*xh;
            yp=scalpr*yh;
            zp=scalpr*zh;
            xp23=xp+x23;
            yp23=yp+y23;
            zp23=zp+z23;
            rp23=(REAL)sqrt((double)(xp23*xp23+yp23*yp23+zp23*zp23));
            xv=xp23/rp23;
            yv=yp23/rp23;
            zv=zp23/rp23;
            x4=x2+BondLen*(sina*xv-cosa*xh);
            y4=y2+BondLen*(sina*yv-cosa*yh);
            z4=z2+BondLen*(sina*zv-cosa*zh);
         }

         hlist->resnum=sNo;
         hlist->insert[0]=sIns;           /* 26.01.96                   */
         padterm(hlist->insert,4);
         strcpy(hlist->atnam,sGHName[4]);
         SetRawAtnam(hlist->atnam_raw, sGHName[4]);   /* 05.12.02       */
         hlist->altpos = ' ';                         /* 03.06.05       */
         hlist->x=x4;
         hlist->y=y4;
         hlist->z=z4;
         ALLOCNEXT(hlist,PDB);
         if(hlist == NULL)
         {
            FREELIST(hlist_start, PDB);
            return(NULL);
         }
         CLEAR_PDB(hlist);

#ifdef DEBUG
         printf("makeh() Type 4 allocated hlist = %d\n",(int)hlist);
#endif

         sNType4++;
      }  /* End of switch                                               */
   }  /* End of HTYPE 1 else clause                                    */

   return(hlist_start);
}

/************************************************************************/
/*>static BOOL AddH(PDB *hlist, PDB **position, int HType)
   -------------------------------------------------------
   Input:   PDB  *hlist       Linked list of hydrogens to be merged
            PDB  **position   Array of PDB pointers for atoms in this
                              residue
            int  HType        Hydrogen type
   Returns: BOOL              Success?

   AddH() merges a list of hydrogens for this atom into the main pdb 
   structure list. Returns FALSE if the procedure failed.

   16.05.90 Original    By: ACRM
   05.12.02 Added setting of atnam_raw
   27.03.03 Fixed memory leak - free the hlist when finished
   03.06.05 Added setting of altpos
*/
static BOOL AddH(PDB *hlist, PDB **position, int HType)
{
   PDB *p,*q,*r,*s;
   int atomcount=0,
       k;

   /* Step through each atom in position list until we find the
      one corresponding to this PGP
   */
   for(k=1;k<16;k++)
   {
      q=hlist;

      if(!position[k]) continue;
      p = position[k];

      /* For PGP types 3 & 5, look for atom in column 3                 */
      if((((HType==3)||(HType==5))
        &&(!(strncmp(p->atnam,sGHName[3],4))))||
        /* For PGP types 1,2,4 look in column 2                         */
        (((HType==1)||(HType==2)||(HType==4))
        &&(!strncmp(p->atnam,sGHName[2],4))))
      {

         /* Copy the atoms from hlist into the PDB list                 */
         s=p;
         r=p->next;           /* Store the pointer to the next record   */
         ALLOCNEXT(p,PDB);    /* Insert a record in the main list       */
         if(p == NULL)
         {
            FREELIST(hlist, PDB);   /* 27.03.03 Fixed memory leak       */
            return(FALSE);
         }
         
         
         p->next=r;           /* Update its pointer                     */
         strcpy(p->record_type,s->record_type);   /* Copy the info into 
                                                     this record        */
         p->atnum = ++atomcount;
         strcpy(p->atnam,q->atnam);
         strcpy(p->atnam_raw, q->atnam_raw);      /* 05.12.02           */
         p->altpos = q->altpos;                   /* 03.06.05           */
         strcpy(p->resnam,s->resnam);
         strcpy(p->chain,s->chain);
         p->resnum=s->resnum;
         strcpy(p->insert,s->insert);
         p->x=q->x;
         p->y=q->y;
         p->z=q->z;
         p->occ=1.0;
         p->bval=20.0;

         if((HType==2)||(HType==3))
         {
            /* For these types there are 2 or three H's to add, so
               add another one
            */
            NEXT(q);
            s=p;
            r=p->next;
            ALLOCNEXT(p,PDB);
            if(p==NULL)
            {
               FREELIST(hlist, PDB);   /* 27.03.03 Fixed memory leak    */
               return(FALSE);
            }
            
            p->next=r;
            strcpy(p->record_type, s->record_type);
            p->atnum = ++atomcount;
            strcpy(p->atnam,q->atnam);
            strcpy(p->atnam_raw, q->atnam_raw);    /* 05.12.02          */
            p->altpos = q->altpos;                 /* 03.06.05          */
            strcpy(p->resnam,s->resnam);
            strcpy(p->chain,s->chain);
            p->resnum=s->resnum;
            strcpy(p->insert,s->insert);
            p->x=q->x;
            p->y=q->y;
            p->z=q->z;
            p->occ=1.0;
            p->bval=20.0;
         }
         if(HType==3)
         {
            /* For this type there are 3 H's to add, so
               add another one
            */
            NEXT(q);
            s=p;
            r=p->next;
            ALLOCNEXT(p,PDB);
            if(p==NULL)
            {
               FREELIST(hlist, PDB);    /* 27.03.03 Fixed memory leak   */
               return(FALSE);
            }
               
            p->next=r;
            strcpy(p->record_type, s->record_type);
            p->atnum = ++atomcount;
            strcpy(p->atnam,q->atnam);
            strcpy(p->atnam_raw, q->atnam_raw);    /* 05.12.02          */
            p->altpos = q->altpos;                 /* 03.06.05          */
            strcpy(p->resnam,s->resnam);
            strcpy(p->chain,s->chain);
            p->resnum=s->resnum;
            strcpy(p->insert,s->insert);
            p->x=q->x;
            p->y=q->y;
            p->z=q->z;
            p->occ=1.0;
            p->bval=20.0;
         }
      }   /* End of matches                                             */
   }  /* End of main list                                               */

   FREELIST(hlist, PDB);    /* 27.03.03 Fixed memory leak               */
   return(TRUE);
}

/************************************************************************/
/*>FILE *OpenPGPFile(char *pgpfile, BOOL AllHyd)
   ---------------------------------------------
   Input:   char   *pgpfile       Name of a PGP file or NULL
            BOOL   AllHyd         If name of PGP not specified, this
                                  flag specified whether all or explicit
                                  hydrogen file required
   Returns: FILE   *              File pointer

   23.08.94 Original    By: ACRM
   28.07.05 Added conditionals for msdos and Mac OS/X
*/
FILE *OpenPGPFile(char *pgpfile, BOOL AllHyd)
{
   char *datadir,
        buffer[160],
        basename[160];
   FILE *fp;
   
   /* If a filename has been specified, just open it and return         */
   if(pgpfile != NULL && pgpfile[0])
   {
      fp = fopen(pgpfile,"r");
      
      return(fp);
   }

   /* No filename specified, so build our own name.
      Select the appropriate base name
   */
   if(AllHyd)
      strcpy(basename,ALLHPGP);
   else
      strcpy(basename,EXPLPGP);
   
   /*** FIXME: This should be changed to use OpenFile() instead       ***/

   /* Try to open file in current directory                             */
   if((fp = fopen(basename,"r")) == NULL)
   {
      /* Failed so build alternative directory/filename                 */
#if (unix || __unix__ || msdos || __msdos__ || __unix || __MACH__ || __APPLE__)
      datadir = getenv(DATAENV);
      
      if(datadir != NULL)
      {
         sprintf(buffer,"%s/%s",datadir,basename);
         fp = fopen(buffer,"r");
      }
#ifdef SCREEN_INFO
      else
      {
         sprintf(buffer,"Error: The %s environment variable has not \
been set.\n",DATAENV);
         screen(buffer);
      }
#endif

#else
      sprintf(buffer,"%s%s",DATADIR,basename);
      fp = fopen(buffer,"r");
#endif
   }
   
   return(fp);
}


/************************************************************************/
/*>static void SetRawAtnam(char *out, char *in)
   --------------------------------------------
   Input:     char    *in    Input string
   Output:    chat    *out   Output string

   Copies the input atom name (a 4-character name) and outputs it to
   a 4-character name starting with a space if the input wtring had
   training spaces. If the input string had no trailing spaces, the
   last character is moved to the front of the string. This creates
   hydrogen atom names in 'raw' format

   05.12.02 Original   By: ACRM
*/
static void SetRawAtnam(char *out, char *in)
{
   char instr[16];
   
   strncpy(instr, in, 15);
   TERMAT(instr, ' ');
   
   if(strlen(instr) > 3)
   {
      out[0] = instr[3];
      out[1] = instr[0];
      out[2] = instr[1];
      out[3] = instr[2];
   }
   else
   {
      out[0]  = ' ';
      strncpy(out+1, instr, 4);
      PADMINTERM(out, 4);
   }
   out[4] = '\0';
}

/************************************************************************/
/*>static PDB *StripDummyH(PDB *pdb, int *nhyd)
   --------------------------------------------
   Input:   PDB    *pdb      Start of PDB linked list
   I/O:     int    *nhyd     Number of hydrogens
   Returns: PDB    *         New PDB linked list

   Strips any dummy hydrogens

   28.11.05 Original   By: ACRM
*/
static PDB *StripDummyH(PDB *pdb, int *nhyd)
{
   PDB *p;

   for(p=pdb; p!=NULL;)
   {
      if((p->atnam[0] == 'H') && 
         (p->x > 9998.0)      &&
         (p->y > 9998.0)      &&
         (p->z > 9998.0))
      {
         DELETE(pdb, p, PDB);
         (*nhyd)--;
      }
      else
      {
         NEXT(p);
      }
   }
   return(pdb);
}

