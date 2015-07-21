/*************************************************************************

   Program:    
   File:       RdSecPDB.c
   
   Version:    V1.3R
   Date:       18.08.98
   Function:   Read secondary structure information from the HELIX, TURN
               and SHEET records in a PDB file
   
   Copyright:  (c) SciTech Software 1990-8
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

**************************************************************************

   Usage:
   ======
   sec = ReadSecPDB(fp,&nsec)

   This routine reads the secondary structure informaton from the
   header of a .PDB file. It reads from file fp and creates a linked
   list of type SECSTRUC. The user should include "pdb.h" 

   On return, the variable nsec, or the value returned by the routine
   should be checked to ensure some secondary structure was returned.
   The secondary structure is stored together with ranges of residues
   to which it applies. The form of the structure should be checked
   in pdb.h

   Input:    fp    *FILE      Pointer to PDB file
   Output:   nsec  *int       Number of sec struc regions identified
   Returns:  sec   *SECSTRUC  Linked list of type SECSTRUC

**************************************************************************

   Revision History:
   =================
   V1.0  22.03.90 Original
   V1.1  20.06.90 Would crash if no secondary structure found---fixed.
   V1.2  09.07.93 Changed allocation scheme. Returns pointer to SECSTRUC.
                  No need to initialise this before calling
   V1.3  18.08.98 Changed SEC to SECSTRUC 'cos of conflict in SunOS

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "SysDefs.h"
#include "MathType.h"

#include "pdb.h"
#include "fsscanf.h"
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
/*>SECSTRUC *ReadSecPDB(FILE *fp, int *nsec)
   -----------------------------------------
   Input:    fp    *FILE      Pointer to PDB file
   Output:   nsec  *int       Number of sec struc regions identified
   Returns:  sec   *SECSTRUC  Linked list of type SECSTRUC

   Reads secondary structure information from the header of a PDB file.
   Returns a pointer to a linked list of type SECSTRUC.

   On return, the variable nsec, or the value returned by the routine
   should be checked to ensure some secondary structure was returned.
   The secondary structure is stored together with ranges of residues
   to which it applies. The form of the structure should be checked
   in pdb.h

   22.03.90 Original
   20.06.90 Would crash if no secondary structure found---fixed.
   09.07.93 Changed allocation scheme. Returns pointer to SECSTRUC.
            No need to initialise this before calling. Uses fsscanf()
   18.08.98 Changed SEC to SECSTRUC 'cos of conflict in SunOS
*/
SECSTRUC *ReadSecPDB(FILE *fp, int *nsec)
{
   int         class,
               res1,
               res2;
   char        buffer[160],
               chain1[8],
               chain2[8],
               ins1[8],
               ins2[8],
               type;
   static char classtab[16] = "HHGGGHGGGH"; /* Class 0 shouldn't occur  */
   SECSTRUC    *p,
               *sec = NULL;
   
   *nsec =  0;

   while(fgets(buffer,159,fp))
   {
      /* Use this as a flag for having found some secondary structure   */
      type    = '\0';
        
      /* Break out of the while loop when we find an ATOM               */
      if(!strncmp(buffer,"ATOM  ",6)) break;
        
      /* Process SHEET records                                          */
      if(!strncmp(buffer,"SHEET ",6))
      {
         fsscanf(buffer,"%21x%1s%4d%1s%5x%1s%4d%1s",
                 chain1, &res1, ins1,
                 chain2, &res2, ins2);
                  
         type = 'E';
      }
        
      /* Process HELIX records                                          */
      if(!strncmp(buffer,"HELIX ",6))
      {
         fsscanf(buffer,"%19x%1s%1x%4d%1s%5x%1s%1x%4d%1s%2d",
                 chain1, &res1, ins1,
                 chain2, &res2, ins2,
                 &class);

         type=classtab[class];
      }
        
      /* Process TURN records                                           */
      if(!strncmp(buffer,"TURN  ",6))
      {
         fsscanf(buffer,"%19x%1s%4d%1s%5x%1s%4d%1s",
                 chain1, &res1, ins1,
                 chain2, &res2, ins2);

         type = 'T';
      }
        
      if(type)    /* We've got some secondary structure                 */
      {
         /* Allcoate space                                              */
         if(sec == NULL)
         {
            INIT(sec,SECSTRUC);
            p = sec;
         }
         else
         {
            ALLOCNEXT(p,SECSTRUC);
         }
         
         /* Check allocation; free list and return if failed            */
         if(p==NULL)
         {
            if(sec != NULL) FREELIST(sec,SECSTRUC);
            *nsec = 0;
            return(NULL);
         }
         
         /* Copy data into linked list                                  */
         strcpy(p->chain1, chain1);
         strcpy(p->chain2, chain2);
         strcpy(p->ins1,   ins1);
         strcpy(p->ins2,   ins2);
         p->res1 = res1;
         p->res2 = res2;
         p->type = type;
         
         (*nsec)++;
      }
   }
   
   /* Return pointer to start of linked list                            */
   return(sec);
}

