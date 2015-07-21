/*************************************************************************

   Program:    
   File:       ErrStack.c
   
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
   StoreError(char *routine, char *error)
   --------------------------------------
   The routine is called with the name of the routine at fault and the
   description of the fault.

   ShowErrors(void *PrintRoutine, BOOL Trace)
   ------------------------------------------
   The routine is called with a pointer to the routine which is to
   do the actual error display and a flag to indicate whether the
   faulty routine names should be displyed. This is only of use if
   the user has access to the source code so should be used for
   debugging purposes only.
   If PrintRoutine is supplied as NULL, the simple PrintAnError()
   routine will be used which displays the error on stderr. More
   complex routines could, for example, show the error in a requester
   or output to a window.
   If a print routine is specified, the routine is called with:
   ShowErrors(void *)MyRoutine, TRUE);

**************************************************************************

   Revision History:
   =================
   V1.0  31.08.94 Original    By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "SysDefs.h"
#include "macros.h"

#include "ErrStack.h"

/************************************************************************/
/* Defines and macros
*/
typedef struct _errorstack
{
   struct _errorstack *next;
   char               *error,
                      *routine;
}  ERRORSTACK;

/************************************************************************/
/* Globals
*/
static ERRORSTACK *sErrorStack = NULL;

/************************************************************************/
/* Prototypes
*/
static void PrintAnError(char *error);

/************************************************************************/
/*>void StoreError(char *routine, char *error)
   -------------------------------------------
   Input:   char *routine        Name of the routine generating the error
            char *error          Description of the error

   Stores an error on the error stack.

   31.08.94 Original    By: ACRM
*/
void StoreError(char *routine, char *error)
{
   static ERRORSTACK *p = NULL;
   
   if(sErrorStack == NULL)
   {
      INIT(sErrorStack, ERRORSTACK);
      p = sErrorStack;
   }
   else
   {
      ALLOCNEXT(p, ERRORSTACK);
   }
   
   if(p != NULL)
   {
      if((p->error = (char *)malloc((strlen(error)+1) * sizeof(char)))
         != NULL)
      {
         strcpy(p->error, error);
      }
      if((p->routine = (char *)malloc((strlen(routine)+1) * sizeof(char)))
         != NULL)
      {
         strcpy(p->routine, routine);
      }
   }
}

/************************************************************************/
/*>void ShowErrors(void *PrintRoutine(char *), BOOL Trace)
   -------------------------------------------------------
   Input:   void  *PrintRoutine(char *)   The print routine or NULL
            BOOL  Trace                   Flag to print routine names

   Display the error stack using the supplied print routine or the
   simple default one if NULL is given.

   31.08.94 Original    By: ACRM
   06.09.94 No longer tries to set PrintRoutine if was NULL (strict ANSI
            compliance)
*/
void ShowErrors(void *PrintRoutine(char *), BOOL Trace)
{
   ERRORSTACK *p;
   char       buffer[160];

   for(p=sErrorStack; p!=NULL; NEXT(p))
   {
      if(Trace)
         sprintf(buffer,"%s :  %s",p->routine,p->error);
      else
         strcpy(buffer,p->error);
      
      if(PrintRoutine == NULL)
         PrintAnError(buffer);
      else
         (*PrintRoutine)(buffer);
   }
}

/************************************************************************/
/*>static void PrintAnError(char *string)
   ---------------------------------------
   Input:   char  *string        A string to be printed

   A simple error printing routine used if NULL given as a parameter
   to ShowErrors()

   31.08.94 Original    By: ACRM
*/
static void PrintAnError(char *string)
{
   fprintf(stderr,"%s\n",string);
}

