/*************************************************************************

   Program:    
   File:       safemem.c
   
   Version:    V1.2
   Date:       03.07.06
   Function:   Safe malloc()/free() routines which check for array 
               overflow on free.
   
   Copyright:  (c) SciTech Software 1995-2006
   Author:     Dr. Andrew C. R. Martin
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   safemalloc() and safefree() are provided as temporary debugging
   replacements for malloc() and free(). They maintain their own
   linked list of malloc()'d memory and allocate a `protection buffer'
   each side of the requested amount of memory. This is filled with a
   given pattern and, when safefree() is called, this buffer is
   checked to ensure the pattern is still present.

   Unlike free(), safefree() is of type BOOL, returning TRUE if
   an error occured.

   SM_SIZE Environment Variable
   The size of the buffer (default 256bytes each side) may be
   controlled by the environment variable SM_SIZE. If code continues
   to core dump, try increasing the value of SM_SIZE.

   SM_FILL Environment Variable 
   The byte used for filling the protection buffers (default 255) is
   defined by the environment variable SM_FILL. If code which core
   dumps with normal malloc()/free(), but not with safemalloc() /
   safefree(), yet no error message is generated, the character which
   is corrupting memory is probably the same as the fill character.
   Try modifying SM_FILL.

**************************************************************************

   Usage:
   ======
   NOTE: safemem.h must be included *after* macros.h

**************************************************************************

   Revision History:
   =================
   V1.0  23.06.95 Original
   V1.1  27.02.98 Added cast to ptr
   V1.2  03.07.06 Added 'ok' to the MEMLIST and added safeleaks()

*************************************************************************/
/* Includes - Note we must *not* include safemem.h since we require the
   standard versions of INITPREV, ALLOCNEXTPREV, etc.
*/
#include <stdio.h>
#include <stdlib.h>
#include "macros.h"
#include "SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
#define CHECKSIZE 256            /* Default protection buffer size      */
#define BLANK     (char)(0xFF)   /* Default blanking character          */
#define MAXBYTE   8              /* Max bytes occupied by an address    */

typedef struct _memlist
{
   char            *start;
   struct _memlist *next,
                   *prev;
   int             length;
   BOOL            ok;
}  MEMLIST;

/************************************************************************/
/* Globals
*/
static MEMLIST *sSafeMemList = NULL;
static int     sCheckSize    = CHECKSIZE;
static char    sBlank        = BLANK;

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>void *safemalloc(int nbytes)
   ----------------------------
   Debugging version of malloc() which creates protection buffers each
   side of the requested memory block.

   23.06.95 Original    By: ACRM
*/
void *safemalloc(int nbytes)
{
   MEMLIST        *p;
   static BOOL    FirstCall = TRUE;
   MEMLIST        *q;
   char           *start;
   int            i;

   /* On first call, check for environment variables                    */
   if(FirstCall)
   {
      char *envvar;
      int  PadValue;
      
      /* Get the protection buffer size (bytes each side of the requested
         chunk of memory
      */
      if((envvar = getenv("SM_SIZE"))!=NULL)
      {
         if(sscanf(envvar,"%d",&sCheckSize))
         {
            /* Ensure sCheckSize is a multiple of MAXBYTE               */
            if(sCheckSize%MAXBYTE)
               sCheckSize = (1 + (int)(sCheckSize/MAXBYTE)) * MAXBYTE;
         }
         else
         {
            fprintf(stderr,"Env. variable, SM_SIZE, is invalid.\n");
         }
      }

      /* Get the padding character to use                               */
      if((envvar = getenv("SM_FILL"))!=NULL)
      {
         if(sscanf(envvar,"%d",&PadValue))
         {
            if(PadValue > 255)
               fprintf(stderr,"Env. variable, SM_FILL, has been \
shortened to the LSB.\n");
            sBlank = PadValue%256;
         }
         else
         {
            fprintf(stderr,"Env. variable, SM_FILL, is invalid.\n");
         }
      }
      
      FirstCall = FALSE;
   }
   
   /* Allocate space in our linked list                                 */
   if(sSafeMemList==NULL)
   {
      INITPREV(sSafeMemList, MEMLIST);
      p=sSafeMemList;
   }
   else
   {
      p=sSafeMemList;
      LAST(p);
      ALLOCNEXTPREV(p, MEMLIST);
   }
   if(p==NULL)
      return(NULL);
   
   /* Allocate the requested memory plus the protection space           */
   if((start=(char *)malloc(nbytes+(2*sCheckSize)))==NULL)
   {
      /* Remove position from our linked list if allocation failed      */
      q = p->prev;
      p->prev->next = NULL;
      free(p);
      p = q;
      return(NULL);
   }
   p->start  = start+sCheckSize;
   p->length = nbytes;
   p->ok     = TRUE;
   
   /* Blank the protection buffer space                                 */
   for(i=0; i<sCheckSize; i++)
   {
      start[i] = sBlank;
      start[i+nbytes+sCheckSize]=sBlank;
   }

   return((void *)(p->start));
}

/************************************************************************/
/*>void safefree(void *ptr)
   ------------------------
   Debugging version of free() which checks protection buffers each
   side of the requested memory block.

   23.06.95 Original    By: ACRM
   27.02.98 Added cast to ptr
*/
BOOL safefree(void *ptr)
{
   MEMLIST *p;
   int     i,
           count=0;

   /* Search our linked list for this address                           */
   for(p=sSafeMemList; p!=NULL; NEXT(p))
   {
      count++;
      
      if(p->start==(char *)ptr)  /* Found it...                         */
      {
         /* Check for array underflow                                   */
         for(i=sCheckSize; i>0; i--)
         {
            if((*(p->start-i) != sBlank))
            {
               fprintf(stderr,"safefree(): Array underflow at \
%lu by %d bytes\n", (ULONG)ptr,i);
               p->ok = FALSE;
               return(TRUE);
            }
         }
         /* Check for array overflow                                    */
         for(i=sCheckSize-1; i>=0; i--)
         {
            if(*(p->start+p->length+i) != sBlank)
            {
               fprintf(stderr,"safefree(): Array overflow at \
%lu by %d bytes\n", (ULONG)ptr,i+1);
               p->ok = FALSE;
               return(TRUE);
            }
         }

         /* No under/over-flow, free the memory and unlink from list    */
         free(p->start-sCheckSize);

         /* Unlink this item from the memory list                       */
         if(p == sSafeMemList)
         {
            NEXT(sSafeMemList);
            if(sSafeMemList)
               sSafeMemList->prev = NULL;
         }
         else
         {
            if(p->prev)
               p->prev->next=p->next;
            if(p->next)
               p->next->prev=p->prev;
         }
         free(p);
         return(FALSE);
      }
   }

   /* Specified address not found in our linked list                    */
   fprintf(stderr,"safefree(): Bad address passed to free: %lu\n",
           (ULONG)ptr);
   fprintf(stderr,"            %d items in memory list\n",count);
   
   return(TRUE);
}

/************************************************************************/
/*>void safeleaks(void)
   --------------------
   Prints a list of any safemalloc()'d memory which was not freed

   03.07.06  Original   By: ACRM
*/
void safeleaks(void)
{
   MEMLIST *p;
   for(p=sSafeMemList; p!=NULL; NEXT(p))
   {
      if(p->ok)
      {
         fprintf(stderr,"safeleaks(): Leaked memory at: %lu of size: \
%lu\n",
                 (long unsigned)p->start, (long unsigned)p->length);
      }
   }
}
