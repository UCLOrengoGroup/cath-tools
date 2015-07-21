/*************************************************************************

   Program:    
   File:       ps.h
   
   Version:    V1.11R
   Date:       23.06.94
   Function:   Include file for PostScript routine
   
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
   Include file for using with PSRoutines.c
   Variables are defined only from the main program. Otherwise they
   are referenced as external.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  06.02.91 Original
   V1.7  25.02.91 Fixed prototypes and definition of PSFile.
   V1.10 07.05.92 Changed all prototypes to doubles
   V1.11 23.06.94 Made gPSFile a global

*************************************************************************/
#ifndef _PS_H
#define _PS_H

/************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "SysDefs.h"
#include "MathType.h"

/************************************************************************/
/* Globals
*/
#ifdef _PS_MAIN   /* ------------------ _PS_MAIN ---------------------- */
   REAL   PSxpicsize    = 5.0,
          PSypicsize    = 5.0,
          PSxoffset     = 1.0,
          PSyoffset     = 2.0;
   FILE   *gPSFile      = NULL;
#else             /* ----------------- Not _PS_MAIN ------------------- */
   extern REAL    PSxpicsize,
                  PSypicsize,
                  PSxoffset,
                  PSyoffset;
   extern FILE    *gPSFile;
#endif            /* -------------------------------------------------- */

/************************************************************************/
/* Prototypes
*/
BOOL PSInit(char *FName, char *creator, char *AltFont);
void PSThick(REAL thickness);
void PSMove(REAL X, REAL Y);
void PSDraw(REAL X, REAL Y);
void PSSetDash(char *linepatt);
void PSClearDash(void);
void PSStroke(void);
void PSFont(char *fontname, REAL size);
void PSLText(REAL X, REAL Y, char *label);
void PSCBText(REAL X, REAL Y, REAL Offset, char *label);
void PSROffText(REAL X, REAL Y, REAL offset, char *label);
void PSLCText(REAL X, REAL Y, char *label);
void PSCTText(REAL X, REAL Y, REAL Offset, char *label);
void PSVText(REAL x, REAL y, REAL xoff, char *text, char *font, REAL size,
             char *label, char *lfont, REAL lsize);
void PSShowText(char *text);
void PSEnd(void);
char *PSCorrectCase(char *font);

#endif
