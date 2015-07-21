/*************************************************************************

   Program:    
   File:       aalist.h
   
   Version:    V3.0
   Date:       18.02.09
   Function:   Include file for amino acid linked lists.
   
   Copyright:  (c) UCL / Dr. Andrew C.R. Martin 2006-2009
   Author:     Dr. Andrew C.R. Martin
   EMail:      andrew@bioinf.org.uk
               
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
   V1.0  21.08.06 Original   By: ACRM
   V3.0  06.11.08 Incorporated into ProFit V3 By: CTP
   V3.0  18.02.09 Moved to bioplib. By: CTP


*************************************************************************/
#ifndef _AALIST_H
#define _AALIST_H


/* Includes
*/
#include "general.h"

/************************************************************************/
/* Defines and macros
*/
typedef struct aa
{
   struct aa *next, *prev;
   int    seqnum;
   BOOL   flag;
   char   res;
} AA;


/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
AA *InsertNextResiduesInAAList(AA *a, char res, int nres);
AA *InsertNextResidueInAAList(AA *a, char res);
char *BuildSeqFromAAList(AA *aa);
AA *InsertResidueInAAListAt(AA *aa, char res, int pos);
AA *InsertResiduesInAAListAt(AA *aa, char res, int nres, int pos);
AA *BuildAAList(char *seq);
int FindAAListOffsetByResnum(AA *aa, int resnum);
AA *FindAAListItemByResnum(AA *aa, int resnum);
void SetAAListFlagByResnum(AA *aa, int resnum);
char *BuildFlagSeqFromAAList(AA *aa, char ch);
int GetAAListLen(AA *aa);

#endif
