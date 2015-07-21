/*************************************************************************

   Program:    
   File:       FNam2PDB.c
   
   Version:    V1.1
   Date:       26.07.95
   Function:   Extract a PDB code from a filename
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1995
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 387 7050 X 3284
   EMail:      INTERNET: martin@biochem.ucl.ac.uk
               
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
   V1.0  24.07.95 Original
   V1.1  26.07.95 Fixed a NULL to '\0'

*************************************************************************/
/* Includes
*/
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
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
/*>char *FNam2PDB(char *filename)
   ------------------------------
   Input:   char *filename     A PDB filename containing a PDB code
   Returns: char *             The PDB code (lower case)
                               NULL if memory allocation fails
                               XXXX if filename blank or NULL or if
                               unable to find PDB code
   This routine attempts to convert a filename stem to a PDB code.

   All the following inputs should produce the same output
   of 1fbj:
   $A:[PDB]PDB1FBJ.ENT
   C:\P1FBJ.PDB
   /pdb/p1fbj.pdb
   /pdb/pdb1fbj.pdb
   1fbj.pdb
   1fbjL
   
   The routine first removes characters from the start of the filename up
   to the last : ] / or \. It then searches for the following possible
   patterns (where N is a digit and X is an alphanumeric)

      pdbNXXX
      pNXXX
      NXXX
      XXXX.pdb
      XXXX.ent

   This should cover just about any filename which includes a legal PDB
   code.

   24.07.95 Original    By: ACRM
   26.07.95 Corrected a NULL to '\0'
*/
char *FNam2PDB(char *filename)
{
   int         length,
               pos;
   char        *temp = NULL,
               *p,
               *start;
   static char pdbcode[8];
   
   /* Check real string given                                           */
   if(filename != NULL)
   {
      /* Check length not zero                                          */
      if((length = strlen(filename))!=0)
      {
         /* Make a temporary copy of the string                         */
         if((temp = (char *)malloc((length+1)*sizeof(char)))==NULL)
            return(NULL);
         strcpy(temp,filename);

         /* Downcase the copy                                           */
         LOWER(temp);

         /* Step from the end of string to char after any : ] / or \    */
         start = temp;
         for(pos=length-1; pos>=0; pos--)
         {
            if((temp[pos] == ':') ||
               (temp[pos] == ']') ||
               (temp[pos] == '/') ||
               (temp[pos] == '\\'))
            {
               start = temp+pos+1;
               break;
            }
         }

         /* Search for `pdb' followed by a digit & 3 alphanumerics      */
         p=start;
         for(;;)
         {
            if((p=strstr(p,"pdb"))!=NULL)
            {
               if(isdigit(p[3]) && 
                  isalnum(p[4]) &&
                  isalnum(p[5]) &&
                  isalnum(p[6]))
               {
                  p[7] = '\0';
                  strcpy(pdbcode, p+3);
                  free(temp);
                  return(pdbcode);
               }
               else
               {
                  p++;
               }
            }
            else
            {
               break;
            }
         }

         /* Search for `p' followed by a digit & 3 alphanumerics        */
         p=start;
         for(;;)
         {
            if((p=strstr(p,"p"))!=NULL)
            {
               if(isdigit(p[3]) &&
                  isalnum(p[4]) &&
                  isalnum(p[5]) &&
                  isalnum(p[6]))
               {
                  p[7] = '\0';
                  strcpy(pdbcode, p+3);
                  free(temp);
                  return(pdbcode);
               }
               else
               {
                  p++;
               }
            }
            else
            {
               break;
            }
         }

         /* Search for digit followed by 3 alphanumerics                */
         for(p=start; *p!='\0'; p++)
         {
            if(isdigit(p[0]) && 
               isalnum(p[1]) &&
               isalnum(p[2]) &&
               isalnum(p[3]))
            {
               p[4] = '\0';
               strcpy(pdbcode, p);
               free(temp);
               return(pdbcode);
            }
         }

         /* Search for 4 characters before .pdb                         */
         if((p=strstr(start,".pdb"))!=NULL)
         {
            p -= 4;
            if(p>=start)
            {
               p[4] = '\0';
               strcpy(pdbcode, p);
               free(temp);
               return(pdbcode);
            }
         }
         
         /* Search for 4 characters before .ent                         */
         if((p=strstr(start,".ent"))!=NULL)
         {
            p -= 4;
            if(p>=start)
            {
               p[4] = '\0'; 
               strcpy(pdbcode, p);
               free(temp);
               return(pdbcode);
            }
         }
         
         /* Give up!                                                    */
         strcpy(pdbcode,"XXXX");
         free(temp);
         return(pdbcode);
      }
   }

   strcpy(pdbcode,"XXXX");
   return(pdbcode);
}
