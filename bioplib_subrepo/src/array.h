/*************************************************************************

   Program:    
   File:       array.h
   
   Version:    V1.5R
   Date:       30.05.02
   Function:   Include file for 2D/3D array functions
   
   Copyright:  (c) SciTech Software 1994-2002
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
   V1.4  18.03.94
   V1.5  30.05.02 Added 3D functions

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

#ifndef _ARRAY_H
#define _ARRAY_H

char **Array2D(int size, int dim1, int dim2);
void FreeArray2D(char **array, int dim1, int dim2);

char ***Array3D(int size, int dim1, int dim2, int dim3);
void FreeArray3D(char ***array, int dim1, int dim2, int dim3);

#endif
