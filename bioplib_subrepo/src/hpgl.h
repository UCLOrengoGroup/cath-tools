/*************************************************************************

   Program:    
   File:       hpgl.h
   
   Version:    V1.0R
   Date:       01.03.94
   Function:   Include file for hpgl
   
   Copyright:  (c) SciTech Software 1991-4
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
   V1.0  25.03.91 Original
*************************************************************************/
#ifndef _HPGL_H
#define _HPGL_H

/************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "SysDefs.h"
#include "MathType.h"

/************************************************************************/
/* Prototypes
*/
BOOL HPGLInit(char *filename, char *AltFont, REAL xmargin, REAL ymargin);
void HPGLPen(int num);
void HPGLMove(REAL x, REAL y);
void HPGLDraw(REAL x, REAL y);
void HPGLSetDash(int style);
void HPGLFont(int font, REAL size);
void HPGLLText(REAL x, REAL y, char *string);
void HPGLCBText(REAL x, REAL y, REAL offset, char *text);
void HPGLROffText(REAL x, REAL y, REAL offset, char *text);
void HPGLLCText(REAL x, REAL y, char *text);
void HPGLCTText(REAL x, REAL y, REAL offset, char *text);
void HPGLVText(REAL x, REAL y, REAL xoff, char *text, int TitleFont, 
               REAL TitleSize, char *label, int LabelFont, REAL LabelSize);
void HPGLEnd(void);
void HPGLShowText(char *text, BOOL orientation, int XBase, int YBase);

#endif
