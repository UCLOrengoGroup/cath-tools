/*************************************************************************

   Program:    
   File:       plotting.h
   
   Version:    V1.0R
   Date:       01.03.94
   Function:   Include file for using plotting routines
   
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
#ifndef _PLOTTING_H
#define _PLOTTING_H

/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "SysDefs.h"
#include "MathType.h"

#include "hpgl.h"
#include "ps.h"

/************************************************************************/
/* Defines
*/
#define DEST_SCREEN  0
#define DEST_PS      1
#define DEST_HPGL    2

/************************************************************************/
/* Prototypes
*/
BOOL AMInitPlot(char *filename, char *title, int dest, REAL OutXSize, 
                REAL OutYSize, REAL OutXOff, REAL OutYOff,
                char *AltFont, REAL xmargin, REAL ymargin,
                REAL DataXMin, REAL DataYMin, REAL DataXMax,
                REAL DataYMax);
void AMSetPen(int dest, int pen);
void AMMove(int dest, REAL x, REAL y);
void AMDraw(int dest, REAL x, REAL y);
void AMSetLineStyle(int dest, int style);
void AMEndLine(int dest);
void AMSetFont(int dest, char *PSFontName, REAL FontSize);
void AMText(int dest, REAL x, REAL y, char *text);
void AMCBText(int dest, REAL x, REAL y, char *text);
void AMRText(int dest, REAL x, REAL y, REAL offset, char *text);
void AMLCText(int dest, REAL x, REAL y, char *text);
void AMCTText(int dest, REAL x, REAL y, REAL CTOffset, char *text);
void AMEndPlot(int dest);
int  PS2HPGLFont(char *font);
char *SimplifyText(char *string);

#endif
