/*************************************************************************

   Program:    
   File:       safemem.h
   
   Version:    V1.0
   Date:       23.06.95
   Function:   Safe malloc()/free() routines which check for array 
               overflow on free.
   
   Copyright:  (c) SciTech Software 1995
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   Fax:        +44 (0) 1372 813069
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
   N.B. This header file should be included *after* macros.h

**************************************************************************

   Revision History:
   =================
   V1.0  23.06.95 Original
   V1.1  03.07.06 Added safeleaks() prototype

*************************************************************************/
#ifndef _SAFEMEM_H
#define _SAFEMEM_H

/* Includes                                                             */
#include "bioplib/SysDefs.h"

/* Prototypes                                                           */
void *safemalloc(int nbytes);
BOOL safefree(void *ptr);
void safeleaks(void);

/* Undefine memory macros defined by macros.h                           */
#ifdef _MACROS_H
#undef INIT
#undef INITPREV
#undef ALLOCNEXT
#undef ALLOCNEXTPREV
#undef FREELIST
#undef DELETE
#endif

/* Redefine macros to use safe versions of malloc()/free()              */
#define INIT(x,y) do { x=(y *)safemalloc(sizeof(y));                     \
                    if(x != NULL) x->next = NULL; } while(0)
#define INITPREV(x,y) do { x=(y *)safemalloc(sizeof(y));                 \
                       if(x != NULL) {x->next = NULL; x->prev = NULL;} } \
                      while(0)
#define ALLOCNEXT(x,y) do { (x)->next=(y *)safemalloc(sizeof(y));        \
                         if((x)->next != NULL) { (x)->next->next=NULL; } \
                         NEXT(x); } while(0)
#define ALLOCNEXTPREV(x,y) do { (x)->next=(y *)safemalloc(sizeof(y));    \
                             if((x)->next != NULL)                       \
                             { (x)->next->prev = (x);                    \
                               (x)->next->next=NULL; }                   \
                               NEXT(x);} while(0)
/* FREELIST takes 2 parameters:
   y: name of list
   z: type of list
*/
#define FREELIST(y,z)   while((y)!=NULL)                                 \
                        {  z *_freelist_macro_q;                         \
                           _freelist_macro_q = (y)->next;                \
                           safefree((char *)(y));                        \
                           (y) = _freelist_macro_q;                      \
                        }
#define ORDFREELIST(y,z)   while((y)!=NULL)                              \
                        {  z *_freelist_macro_q;                         \
                           _freelist_macro_q = (y)->next;                \
                           free((char *)(y));                            \
                           (y) = _freelist_macro_q;                      \
                        }

/*>DELETE(start, item, type)
   -------------------------
   Deletes item from a linked list
   start may be modified if item is the first in the list
   item is returned as the pointer to the next item in the list (i.e.
   as item->next). One can therefore simply call the routine N time
   to delete N items. If start or item is NULL, does nothing

   16.02.95 Original    By: ACRM
*/
#define DELETE(x, y, z)                                                  \
do {                                                                     \
   z *_delete_macro_p,                                                   \
     *_delete_macro_prev = NULL,                                         \
     *_delete_macro_temp,                                                \
     *_delete_macro_temp2;                                               \
   if((x)!=NULL && (y)!=NULL)                                            \
   {                                                                     \
      for(_delete_macro_p=(x);                                           \
          _delete_macro_p!=NULL;                                         \
          NEXT(_delete_macro_p))                                         \
      {                                                                  \
         if(_delete_macro_p == (y))                                      \
         {                                                               \
            _delete_macro_temp2 = (y)->next;                             \
            if(_delete_macro_prev == NULL)                               \
            {                                                            \
               _delete_macro_temp = (x)->next;                           \
               safefree(x);                                              \
               (x) = _delete_macro_temp;                                 \
               break;                                                    \
            }                                                            \
            else                                                         \
            {                                                            \
               _delete_macro_prev->next = _delete_macro_p->next;         \
               safefree(_delete_macro_p);                                \
            }                                                            \
         }                                                               \
         _delete_macro_prev = _delete_macro_p;                           \
      }                                                                  \
      (y) = _delete_macro_temp2;                                         \
   }                                                                     \
}  while(FALSE)

#define SAFEINIT(x,y) INIT(x,y)
#define SAFEINITPREV(x,y) INITPREV(x,y)
#define SAFEALLOCNEXT(x,y) ALLOCNEXT(x,y)
#define SAFEALLOCNEXTPREV(x,y) ALLOCNEXTPREV(x,y)
#define SAFEFREELIST(y,z) FREELIST(y,z)
#define SAFEDELETE(x, y, z) DELETE(x, y, z)

#endif
