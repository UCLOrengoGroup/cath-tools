/*************************************************************************

   Program:    
   File:       TrueSeqLen.c
   
   Version:    V1.10
   Date:       02.10.00
   Function:   
   
   Copyright:  (c) SciTech Software 1993-2000
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
   V1.0  29.09.92 Original
   V1.1  07.06.93 Corrected allocation
   V1.2  18.06.93 Handles multi-chains and skips NTER and CTER residues.
                  Added SplitSeq()
   V1.3  09.07.93 SplitSeq() cleans up properly if allocation failed
   V1.4  11.05.94 Added TrueSeqLen()
   V1.5  13.05.94 Fixed bug in PDB2Seq().
                  Added KnownSeqLen().
   V1.6  07.09.94 Fixed allocation bug in SplitSeq()
   V1.7  19.07.95 Added check for ATOM records
   V1.8  24.01.96 Fixed bug when no ATOM records in linked list
                  Returns a blank string
   V1.9  26.08.97 Renamed DoPDB2Seq() with handling of Asx/Glx and
                  protein-only. Added macros to recreate the
                  old PDB2Seq() interface and similar new calls
   V1.10 02.10.00 Added NoX option

*************************************************************************/
/* Includes
*/

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
/*>int TrueSeqLen(char *sequence)
   ------------------------------
   Input:   char  *sequence    A sequence containing deletions
   Returns: int                Length without deletions

   Scans a 1-letter code sequence and calculate the length without
   `-' or ` ' residues

   14.04.94 Original    By: ACRM
*/
int TrueSeqLen(char *sequence)
{
   int length = 0,
       i = 0;
   
   for(i=0; sequence[i]; i++)
   {
      if(sequence[i] != '-' && sequence[i] != ' ')
         length++;
   }
   
   return(length);
}

