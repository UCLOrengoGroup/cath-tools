/*************************************************************************

   Program:    
   File:       GetCrystPDB.c
   
   Version:    V1.0
   Date:       12.10.95
   Function:   
   
   Copyright:  (c) SciTech Software 1993-1995
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
   V1.0R 12.10.05 Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include "SysDefs.h"
#include "MathType.h"
#include "pdb.h"
#include "fsscanf.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/


/************************************************************************/
/*>int GetCrystPDB(FILE *fp, VEC3F *UnitCell, VEC3F *CellAngles,
                   char *spacegroup,
                   REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4])
   --------------------------------------------------------------
   Input:   FILE    *fp                 Input file pointer
   Output:  VEC3F   *UnitCell           The unit cell dimensions
            VEC3F   *CellAngles         The unit cell angles
            char    *spacegroup         The spacegroup
            REAL    OrigMatrix[3][4]    The origin matrix
            REAL    ScaleMatrix[3][4]   The scale matrix
   Returns: int                         Flags for elements read:
                                        0: Nothing at all
                                        XTAL_DATA_CRYST: Unit cell
                                        XTAL_DATA_ORIGX: Origin matrix
                                        XTAL_DATA_SCALE: Scale matrix
                                        (These are ORed together)

   Read the crystal parameters (unit cell, spacegroup, origin and scale 
   matrices) out of a PDB file.

   Stops searching as soon as an ATOM or HETATM record is hit and leaves
   the file in a state ready for ReadPDB() to do it's stuff (i.e. with
   the current file pointer at the first ATOM or HETATM record).

   12.10.95 Original    By: ACRM
   17.10.95 Correctly returns retval if no ATOM/HETATM records not found
*/
int GetCrystPDB(FILE *fp, VEC3F *UnitCell, VEC3F *CellAngles,
                char *spacegroup,
                REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4])
{
   int  ch,
        i, j,
        record,
        retval = 0x0000;
   char buffer[MAXBUFF];
   BOOL FirstChar = TRUE,
        WholeLine = FALSE;

   /* Initialise matrices and cell dimensions                           */
   UnitCell->x   = UnitCell->y   = UnitCell->z   = (REAL)1.0;
   CellAngles->x = CellAngles->y = CellAngles->z = (REAL)90.0*PI/180.0;
   for(i=0; i<3; i++)
   {
      for(j=0; j<4; j++)
      {
         OrigMatrix[i][j]  = (REAL)0.0;
         ScaleMatrix[i][j] = (REAL)0.0;
      }
      OrigMatrix[i][i] = (REAL)1.0;
   }
   strcpy(spacegroup,"P");

   while(((ch=fgetc(fp)) != EOF) && !feof(fp))
   {
      if(ch == '\n')
      {
         FirstChar = TRUE;
      }
      else if(FirstChar)
      {
         WholeLine = FALSE;
         
         if(ch=='A' || ch=='C' || ch=='H' || ch=='O' || ch=='S')
         {
            buffer[0] = (char)ch;
            for(i=1; i<6; i++)
            {
               ch = fgetc(fp);
               if(ch==EOF || feof(fp))
                  return(FALSE);
                  
               buffer[i] = (char)ch;
            }
            buffer[6] = '\0';
            
            /* See if it was an ATOM or HETATM record                   */
            if(!strncmp(buffer,"ATOM  ",6) || !strncmp(buffer,"HETATM",6))
            {
               /* If so, shove the characters back into the input stream
                  and return
               */
               for(i=5; i>=0; i--)
                  ungetc(buffer[i], fp);
               return(retval);
            }

            /* See if it was a CRYST record                             */
            if(!strncmp(buffer,"CRYST1",6))
            {
               /* If so, shove the characters back into the input stream,
                  extract the crystal parameters and set the return value
               */
               for(i=5; i>=0; i--)
                  ungetc(buffer[i], fp);
               fgets(buffer, MAXBUFF, fp);
               WholeLine = TRUE;
               fsscanf(buffer,
                       "%6x%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf%1x%14s",
                       &(UnitCell->x), 
                       &(UnitCell->y), 
                       &(UnitCell->z),
                       &(CellAngles->x), 
                       &(CellAngles->y), 
                       &(CellAngles->z),
                       spacegroup);
               CellAngles->x *= PI/180.0;
               CellAngles->y *= PI/180.0;
               CellAngles->z *= PI/180.0;
              
               retval |= XTAL_DATA_CRYST;
            }               

            /* See if it was a ORIGX record                             */
            if(!strncmp(buffer,"ORIGX",5))
            {
               /* If so, shove the characters back into the input stream,
                  extract the crystal parameters and set the return value
               */
               for(i=5; i>=0; i--)
                  ungetc(buffer[i], fp);
               fgets(buffer, MAXBUFF, fp);
               fsscanf(buffer,"%5x%1d",&record);
               WholeLine = TRUE;
               record--;
               fsscanf(buffer,
                       "%10x%10.6lf%10.6lf%10.6lf%15.5lf",
                       &(OrigMatrix[record][0]),
                       &(OrigMatrix[record][1]),
                       &(OrigMatrix[record][2]),
                       &(OrigMatrix[record][3]));
              
               retval |= XTAL_DATA_ORIGX;
            }               

            /* See if it was a SCALE record                             */
            if(!strncmp(buffer,"SCALE",5))
            {
               /* If so, shove the characters back into the input stream,
                  extract the crystal parameters and set the return value
               */
               for(i=5; i>=0; i--)
                  ungetc(buffer[i], fp);
               fgets(buffer, MAXBUFF, fp);
               fsscanf(buffer,"%5x%1d",&record);
               WholeLine = TRUE;
               record--;
               fsscanf(buffer,
                       "%10x%10.6lf%10.6lf%10.6lf%15.5lf",
                       &(ScaleMatrix[record][0]),
                       &(ScaleMatrix[record][1]),
                       &(ScaleMatrix[record][2]),
                       &(ScaleMatrix[record][3]));
              
               retval |= XTAL_DATA_SCALE;
            }               

         }
         
         FirstChar = WholeLine?TRUE:FALSE;
      }
   }
   return(retval);
}


