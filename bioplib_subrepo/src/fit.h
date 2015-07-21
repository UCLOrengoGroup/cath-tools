/*************************************************************************

   Program:    
   File:       fit.h
   
   Version:    V1.1R
   Date:       01.03.94
   Function:   Include file for least squares fitting
   
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
   V1.0  04.02.91 Original
   V1.1  08.12.92 Removed qikfit() prototype as is static

*************************************************************************/
#ifndef _FIT_H
#define _FIT_H

#include "MathType.h"
#include "SysDefs.h"

/* Prototypes for functions defined in fit.c                            */
BOOL matfit(COOR *x1, COOR *x2, REAL rm[3][3], int n, REAL *wt1, 
            BOOL column);

#endif

