/*************************************************************************

   Program:    
   File:       RdSSPDB.c
   
   Version:    V1.0R
   Date:       01.03.94
   Function:   Read disulphide information from header records of a PDB
               file
   
   Copyright:  (c) SciTech Software 1993-4
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
   Read the header records of a PDB file to find disulphide information.
   Builds a linked list of type DISULPHIDE.

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
#include <stdlib.h>

#include "SysDefs.h"

#include "pdb.h"
#include "macros.h"
#include "fsscanf.h"

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
/*>DISULPHIDE *ReadDisulphidesPDB(FILE *fp, BOOL *error)
   -----------------------------------------------------
   Input:   FILE        *fp     PDB file pointer
   Output:  BOOL        *error  Success
   Returns: DISULPHIDE  *       Linked list of disulphide information.
                                NULL if none found or error (Check flag)

   Searches a PDB file for SSBOND records and constructs a linked list
   of information from these records.
   Returns NULL if no disulphide information found. If memory allocation
   fails, the DISULPHIDE linked list formed thus far is returned and 
   the error flag is set to TRUE

   14.10.93 Original   By: ACRM
*/
DISULPHIDE *ReadDisulphidesPDB(FILE *fp, BOOL *error)
{
   DISULPHIDE *dis = NULL,
              *p   = NULL;
   char       buffer[160];

   *error = FALSE;

   while(fgets(buffer,160,fp))
   {
      /* Exit as soon as we reach an ATOM record                        */
      if(!strncmp(buffer,"ATOM  ",6)) break;
      if(!strncmp(buffer,"SSBOND",6))
      {
         /* Allocate memory in linked list                              */
         if(dis==NULL)
         {
            INIT(dis,DISULPHIDE);
            p=dis;
         }
         else
         {
            ALLOCNEXT(p,DISULPHIDE);
         }

         if(p==NULL)
         {
            *error = TRUE;
            return(dis);
         }

         /* Read data out of SSBOND record                              */
         fsscanf(buffer,"%15x%1s%5d%1s%7x%1s%5d%1s",
                 p->chain1,&p->res1,p->insert1,
                 p->chain2,&p->res2,p->insert2);
      }
   }

   return(dis);
}

