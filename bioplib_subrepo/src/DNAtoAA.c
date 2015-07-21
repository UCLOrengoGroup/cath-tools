/*************************************************************************

   Program:    
   File:       DNAtoAA.c
   
   Version:    V1.0R
   Date:       11.05.94
   Function:   Convert DNA codons to amino acid 1-letter code
   
   Copyright:  (c) SciTech Software 1994
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
   #include "bioplib/seq.h" to define prototype

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include "macros.h"
#include <string.h>
#include <ctype.h>

/************************************************************************/
/* Defines and macros
*/
#define NUCINDEX(x) ((x)=='T' ? (0) : \
                    ((x)=='U' ? (0) : \
                    ((x)=='C' ? (1) : \
                    ((x)=='A' ? (2) : \
                     (3)))))

/************************************************************************/
/* Globals
*/
static char *sAACode[4][4] =
{  { "FFLL", "SSSS", "YYXX", "CCXW" },   /* TTX, TCX, TAX, TGX          */
   { "LLLL", "PPPP", "HHQQ", "RRRR" },   /* CTX, CCX, CAX, CGX          */
   { "IIIM", "TTTT", "NNKK", "SSRR" },   /* ATX, ACX, AAX, AGX          */
   { "VVVV", "AAAA", "DDEE", "GGGG" }    /* GTX, GCX, GAX, GGX          */
} ;

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>char DNAtoAA(char *dna)
   -----------------------
   Input:   char  *dna        DNA/RNA codon
   Returns: char              1-letter amino acid code (X=termination)

   Converts a nucleic acid codon to the 1-letter amino acid equivalent.
   Termination codons are returned as X. No special action is taken
   for initiation codons.

   18.04.94 Original    By: ACRM
*/
char DNAtoAA(char *dna)
{
   char buffer[8], *p;
   int idx1, idx2, idx3;
   
   KILLLEADSPACES(p,dna);
   
   strncpy(buffer,p,8);
   buffer[7] = '\0';
   UPPER(buffer);

   idx1 = NUCINDEX(buffer[0]);
   idx2 = NUCINDEX(buffer[1]);
   idx3 = NUCINDEX(buffer[2]);
   
   return(sAACode[idx1][idx2][idx3]);
}

