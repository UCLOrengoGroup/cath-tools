/*************************************************************************

   Program:    
   File:       ErrStack.h
   
   Version:    V1.0R
   Date:       31.08.94
   Function:   Build and print an error stack for program failure.
   
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

   This set of routines allows a stack of errors to be created. 
   When a program has a fatal error, the StoreError() routine is called
   to place the error on the stack. As the program un-winds, each
   routine which fails stores it's error. Finally, before the program
   actually exits, it calls ShowErrors() to display the error stack.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  31.08.94 Original    By: ACRM

*************************************************************************/
#ifndef _ERRSTACK_H
#define _ERRSTACK_H
/************************************************************************/
/* Includes
*/
#include "SysDefs.h"

/************************************************************************/
/* Prototypes
*/
void StoreError(char *routine, char *error);
void ShowErrors(void *PrintRoutine(char *), BOOL Trace);

/************************************************************************/
#endif
