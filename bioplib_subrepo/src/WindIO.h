/*************************************************************************

   Program:    
   File:       WindIO.h
   
   Version:    V1.3R
   Date:       18.10.95
   Function:   Header for window/normal interface routines

   Copyright:  (c) SciTech Software 1993-5
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
#ifndef _WINDIO_H
#define _WINDIO_H

#include "SysDefs.h"

void screen(char *string);
void prompt(char *string);
void RePrompt(void);
void GetKybdString(char *string, int maxlen);
void PagingOn(void);
void PagingOff(void);
void WindowMode(BOOL mode);
void WindowInteractive(BOOL mode);
int YorN(char deflt);

#endif
