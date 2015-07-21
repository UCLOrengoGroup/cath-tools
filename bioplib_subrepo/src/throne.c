/*************************************************************************

   Program:    
   File:       throne.c
   
   Version:    V1.7
   Date:       18.02.09
   Function:   Convert between 1 and 3 letter aa codes
   
   Copyright:  (c) SciTech Software 1993-2009
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  29.09.92 Original    By: ACRM   
   V1.1  11.03.94 Added PCA, ASX and GLX to translation table.
                  PCA translates to E
                  Added routines to handle asx/glx
   V1.2  25.07.95 handles nucleic acids
                  Sets the gBioplibSeqNucleicAcid flag if it's a 
                  nucleic acid.
   V1.3  08.03.07 Added PGA (Pyroglutamate) to translation table 
                  (same as PCA: pyrrolidone carboxylic acid).
                  Note that it isn't clear whether this should translate
                  to Glu or Gln
   V1.4  21.07.08 Added CGN (5-OXO-PYRROLIDINE-2-CARBALDEHYDE) which
                  is again the same as PCA
   V1.5  19.12.08 Corrected NUMAAKNOWN - wasn't looking at U or X as
                  these were > NUMAAKNOWN!
   V1.6  04.02.09 onethr() was not properly working from end of list
                  for nucleic acids
   V1.7  18.02.09 Fixed for new PDB files which have "  DT" etc for DNA
                  sequences

*************************************************************************/
/* Includes
*/
#include <string.h>
#include "SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
#define NUMAAKNOWN 37

/************************************************************************/
/* Globals
*/

/* N.B. The order in sTab1[] and sTab3[] must be the same and they must
   end with X/UNK.
   Also, nucleic acids must come *after* amino acids.
*/
/* Don't forget to fix NUMAAKNOWN if adding to this table!              */
static char sTab1[]    = {'A','C','D','E','F',
                          'G','H','I','K','L',
                          'M','N','P','Q','R',
                          'S','T','V','W','Y',
                          'E','B','Z','E','E',
                          'A','T','C','G','U','I',
                          'A','T','C','G','I','X'
                         };
/* Don't forget to fix NUMAAKNOWN if adding to this table!              */
static char sTab3[][8] = {"ALA ","CYS ","ASP ","GLU ","PHE ",
                          "GLY ","HIS ","ILE ","LYS ","LEU ",
                          "MET ","ASN ","PRO ","GLN ","ARG ",
                          "SER ","THR ","VAL ","TRP ","TYR ",
                          "PCA ","ASX ","GLX ","PGA ","CGN ",
                          "  A ","  T ","  C ","  G ","  U ","  I ",
                          " DA "," DT "," DC "," DG "," DI ","UNK "
                         };
/* Don't forget to fix NUMAAKNOWN if adding to this table!              */

BOOL gBioplibSeqNucleicAcid = FALSE;

/************************************************************************/
/* Prototypes
*/


/************************************************************************/
/*>char throne(char *three)
   ------------------------
   Input:   char  *three    Three letter code
   Returns: char            One letter code

   Converts 3-letter code to 1-letter code.
   Handles ASX and GLX as X
   
   29.09.92 Original    By: ACRM
   11.03.94 Modified to handle ASX and GLX in the tables
   25.07.95 Added handling of gBioplibSeqNucleicAcid
*/
char throne(char *three)
{
   int j;

   if(three[0] == ' ' && three[1] == ' ')
      gBioplibSeqNucleicAcid = TRUE;
   else
      gBioplibSeqNucleicAcid = FALSE;

   if(three[2] == 'X')
      return('X');

   for(j=0;j<NUMAAKNOWN;j++)
      if(!strncmp(sTab3[j],three,3)) return(sTab1[j]);

   /* Only get here if the three letter code was not found              */
   return('X');
}


/************************************************************************/
/*>char thronex(char *three)
   -------------------------
   Input:   char  *three    Three letter code
   Returns: char            One letter code

   Converts 3-letter code to 1-letter code.
   Handles ASX and GLX as B and Z.
   
   29.09.92 Original    By: ACRM
   25.07.95 Added handling of gBioplibSeqNucleicAcid
*/
char thronex(char *three)
{
   int j;

   if(three[0] == ' ' && three[1] == ' ')
      gBioplibSeqNucleicAcid = TRUE;
   else
      gBioplibSeqNucleicAcid = FALSE;

   for(j=0;j<NUMAAKNOWN;j++)
      if(!strncmp(sTab3[j],three,3)) return(sTab1[j]);

   /* Only get here if the three letter code was not found              */
   return('X');
}


/************************************************************************/
/*>char *onethr(char one)
   ----------------------
   Input:   char  one     One letter code
   Returns: char  *       Three letter code (padded to 4 chars with a 
                          space)

   Converts 1-letter code to 3-letter code (actually as 4 chars).

   07.06.93 Original    By: ACRM
   25.07.95 If the gBioplibSeqNucleicAcid flag is set, assumes nucleic
            acids rather than amino acids
   03.02.09 Fixed nucleic search - j was incrementing instead of 
            decrementing!
*/
char *onethr(char one)
{
   int j;

   if(gBioplibSeqNucleicAcid) /* Work from end of table                 */
   {
      for(j=NUMAAKNOWN-1;j>=0;j--)
         if(sTab1[j] == one) return(sTab3[j]);
   }
   else                       /* Work from start of table               */
   {
      for(j=0;j<NUMAAKNOWN;j++)
         if(sTab1[j] == one) return(sTab3[j]);
   }

   /* Only get here if the one letter code was not found                */
   return(sTab3[NUMAAKNOWN-1]);
}

