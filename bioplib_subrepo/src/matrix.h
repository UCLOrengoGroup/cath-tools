/*************************************************************************

   Program:    
   File:       matrix.h
   
   Version:    V1.6R
   Date:       27.09.95
   Function:   Include file for matrix operations
   
   Copyright:  (c) SciTech Software 1995
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

#ifndef _MATRIX_H
#define _MATRIX_H
#include "MathType.h"

void MatMult3_33(VEC3F vecin, REAL matin[3][3], VEC3F *vecout);
void MatMult33_33(REAL a[3][3], REAL b[3][3], REAL out[3][3]);
void invert33(REAL s[3][3], REAL ss[3][3]);
void CreateRotMat(char direction, REAL angle, REAL matrix[3][3]);
REAL VecDist(REAL *a, REAL *b, int len);

#endif
