/*************************************************************************

   Program:    
   File:       Macros.h
   
   Version:    V2.17
   Date:       10.04.08
   Function:   Useful macros
   
   Copyright:  SciTech Software 1991-2008
   Author:     Andrew C. R. Martin
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
   If not Amiga defines abs().
   Defines max(), min() and PI if not done.
   Defines list handling macros.
   Defines newline() and toggle() macros.
   
**************************************************************************

   Usage:
   ======
   INIT(x,y)            Initialise list of name x and type y. 
                        Set x->next to NULL
   INITPREV(x,y)        Ditto, but also sets x->prev to NULL
   NEXT(x)              Step on in linked list
   PREV(x)              Step back in linked list
   ALLOCNEXT(x,y)       Allocate next item in list and step on
   ALLOCNEXTPREV(x,y)   Allocate next item in list and step on. 
                        Also set ->prev item in next item
   LAST(x)              Move to end of list
   FREELIST(y,z)        Free list y of type z
   NEWLINE              Print a newline character ty stdout
   TOGGLE(x)            Toggle a flag
   RANGECHECK(x,y,z)    Return x constrained to range y to z
   TERMINATE(x)         Terminate a string at the first \n
   MAX(x,y)             max() as macro
   MIN(x,y)             min() as macro
   ABS(x,y)             abs() as macro
   UPPER(x)             Converts a string to upper case
   KILLLEADSPACES(x,y)  Makes x a pointer into string y after any spaces
                        or tabs.
   D(BUG)               Prints the BUG string if DEBUG is defined first
   DELETE(lst,itm,type) Deletes (itm) from linked list (lst) of type
                        (type)
   TESTINARRAY(x,l,y,r) Tests whether value (y) is in array (x) if length
                        (l) returning the result in (r)
   FINDINARRAY(x,l,y,r) Finds value (y) is in array (x) if length
                        (l) returning the offset in (r). Offset is -1 if
                        not found
   SET(x,y)             Sets bit y (a hex value) in variable x
   UNSET(x,y)           Clears bit y (a hex value) in variable x
   ISSET(x,y)           Tests bit y (a hex value) in variable x
   TERMAT(x,y)          Terminates character string x at first character y
   KILLTRAILSPACES(x)   Terminate string to remove any trailing white 
                        space
   PROMPT(fp,x)         Issue a prompt to stdout if fp is a terminal
   PADMINTERM(str,len)  Pads a string to len chars only if it is shorter
   DELETEDOUBLE(lst,itm,type) Deletes (itm) from a doubly linked list 
                        (lst) of type (type)
   PADCHARMINTERM(str,char,len) Pads a string to len chars with specified
                        character only if it is shorter
   DOTIFY(str)          Replace ' ' with '.' in string
   DEDOTIFY(str)        Replace '.' with ' ' in string
   FINDPREV(p, start, q) Set p to item in linked list start before q
   
**************************************************************************

   Revision History:
   =================
   V1.0  06.02.91 Original
   V1.1  15.02.91 Moved PI definition to non-Amiga's only
   V1.2  21.03.91 Added RANGECHECK
   V1.3  06.09.91 Added DIST, DISTSQ and Vec3f
   V1.4  09.09.91 Fixed multi-command macros with {}
   V1.5  24.01.92 Fixed for 32 bit addresses and added malloc checks.
   V1.6  03.04.92 Small change to ALLOCNEXT and ALLOCNEXTPREV, so
                  will do a NEXT() even if malloc() fails.
   V1.7  06.05.92 Added TERMINATE()
   V1.8  06.07.92 Added MAX(), MIN() and ABS()
   V1.9  22.07.92 Fixed ABS()
   V1.10 28.09.92 Added TRUE & FALSE and UPPER()
   V1.11 03.11.92 Changed TOGGLE and newline is now upper case.   
   V1.12 16.11.92 Added KILLLEADSPACES()
   V1.13 18.11.92 Fixed UPPER() for MicrosoftC which returns strlen()
                  as unsigned
   V1.14 20.11.92 ABS() now uses 0 rather than 0.0, so we don't
                  try to use floats with ints...
   V2.0  24.11.92 Removed all small letter macros
   V2.1  12.07.93 Added double include check, moved math definitions
                  to MathType.h and added LOWER()
   V2.2  07.10.93 UPPER() and LOWER() check case first for ESV
                  compatibility
   V2.3  23.05.94 Added D(BUG)
   V2.4  14.07.94 Added do{}while(0) bracketing of all multi-line macros
   V2.5  21.11.94 ABS, MAX and MIN check that they're not already defined
   V2.6  16.02.95 Added DELETE()
   V2.7  21.02.95 Updated some internal variable names
   V2.8  02.08.95 Added TESTINARRAY(), FINDINARRAY(), 
                  SET(), UNSET() and ISSET()
   V2.9  20.11.95 Added TERMAT()
   V2.10 06.02.96 Added KILLTRAILSPACES()
   V2.11 14.06.96 Added PROMPT()
   V2.12 23.07.96 Added PADMINTERM()
   V2.13 19.09.96 Include ctype for UPPER() etc
   V2.14 13.03.99 Added DELETEDOUBLE()
   V2.15 01.03.01 Added DOTIFY() DEDOTIFY() PADCHARMINTERM() SUBSCHAR()
   V2.16 25.01.06 Added FINDPREV()
   V2.17 10.04.08 Fixed bug in DELETE() - the break was not properly
                  stopping prev from being changed

*************************************************************************/
#ifndef _MACROS_H
#define _MACROS_H

/***************************** Includes *********************************/
#include <ctype.h>

/**************************** Definitions *******************************/
#ifndef PI
#define PI (4.0 * atan(1.0))
#endif

#ifndef TRUE
#define TRUE  1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/***************************** Maths macros *****************************/
#define RANGECHECK(x,y,z) ((x)<(y)) ? (y) : ((x)>(z)) ? (z) : (x)
#define DISTSQ(a,b) (((a)->x - (b)->x) * ((a)->x - (b)->x) + \
                     ((a)->y - (b)->y) * ((a)->y - (b)->y) + \
                     ((a)->z - (b)->z) * ((a)->z - (b)->z))
#define DIST(a,b) sqrt(((a)->x - (b)->x) * ((a)->x - (b)->x) + \
                       ((a)->y - (b)->y) * ((a)->y - (b)->y) + \
                       ((a)->z - (b)->z) * ((a)->z - (b)->z))
#ifndef ABS
#define ABS(x)   (((x)<0)   ? (-(x)) : (x))
#endif

#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a)    : (b))
#define MIN(a,b) (((a)<(b)) ? (a)    : (b))
#endif

/***************************** List macros ******************************/
#define INIT(x,y) do { x=(y *)malloc(sizeof(y)); \
                    if(x != NULL) x->next = NULL; } while(0)
#define INITPREV(x,y) do { x=(y *)malloc(sizeof(y));\
                        if(x != NULL) {x->next=NULL; x->prev=NULL;} } \
                      while(0)
#define NEXT(x) (x)=(x)->next
#define PREV(x) (x)=(x)->prev
#define ALLOCNEXT(x,y) do { (x)->next=(y *)malloc(sizeof(y));\
                         if((x)->next != NULL) { (x)->next->next=NULL; }\
                         NEXT(x); } while(0)
#define ALLOCNEXTPREV(x,y) do { (x)->next=(y *)malloc(sizeof(y));\
                             if((x)->next != NULL)\
                             { (x)->next->prev = (x); \
                               (x)->next->next=NULL; }\
                               NEXT(x);} while(0)
#define LAST(x)   while((x)->next != NULL) NEXT(x)
/* FREELIST takes 2 parameters:
   y: name of list
   z: type of list
*/
#define FREELIST(y,z)   while((y)!=NULL) \
                        {  z *_freelist_macro_q; \
                           _freelist_macro_q = (y)->next; \
                           free((char *)(y)); \
                           (y) = _freelist_macro_q; \
                        }

/*>DELETE(start, item, type)
   -------------------------
   Deletes (item) from a linked list.
   (start) will be modified if (item) is the first in the list.
   (item) is returned as the pointer to the next item in the list (i.e.
   as item->next). One can therefore simply call the routine N times
   to delete N items. If (start) or (item) is NULL, does nothing

   16.02.95 Original    By: ACRM
   10.04.08 Fixed position of break. By: CTP
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
               free(x);                                                  \
               (x) = _delete_macro_temp;                                 \
            }                                                            \
            else                                                         \
            {                                                            \
               _delete_macro_prev->next = _delete_macro_p->next;         \
               free(_delete_macro_p);                                    \
            }                                                            \
            break;                                                       \
         }                                                               \
         _delete_macro_prev = _delete_macro_p;                           \
      }                                                                  \
      (y) = _delete_macro_temp2;                                         \
   }                                                                     \
}  while(FALSE)


/*>DELETEDOUBLE(start, item, type)
   -------------------------------
   Deletes (item) from a doubly linked list.
   (start) will be modified if (item) is the first in the list.
   (item) is returned as the pointer to the next item in the list (i.e.
   as item->next). One can therefore simply call the routine N times
   to delete N items. If (start) or (item) is NULL, does nothing

   13.03.99 Original    By: ACRM
*/
#define DELETEDOUBLE(s, x, y)                                            \
        do { y *_deleteandnext_macro_temp;                               \
             if(((s)!=NULL) && ((x)!=NULL))                              \
             {  if((x)==(s)) (s) = (x)->next;                            \
                _deleteandnext_macro_temp = (x)->next;                   \
                if((x)->prev != NULL) (x)->prev->next = (x)->next;       \
                if((x)->next != NULL) (x)->next->prev = (x)->prev;       \
                free(x);                                                 \
                (x) = _deleteandnext_macro_temp;                         \
        }  } while(0)

/*>FINDPREV(ptr, start, item)
   --------------------------
   Searches a linked list beginning at (start) to find the item which
   preceeds (item). Its address is put into (ptr). If (item) is the
   same as (start) or (item) is not found, then the routine returns
   NULL in (ptr)
   This is used when wanting to look at the previous item in a singly
   linked list.

   26.01.06 Original    By: ACRM
*/
#define FINDPREV(p, s, l)                                                \
        do { p = (s);                                                    \
             if((s)==(l))                                                \
             { p = NULL; } else                                          \
             {                                                           \
               while((p != NULL) && (p->next != (l)))                    \
               {  p = p->next;                                           \
           } } } while(0)


/***************************** Misc. macros *****************************/
#define NEWLINE printf("\n")

#define TOGGLE(x) (x) = (x) ? FALSE : TRUE

#define TERMINATE(x) do {  int _terminate_macro_j;                    \
                        for(_terminate_macro_j=0;                     \
                            (x)[_terminate_macro_j];                  \
                            _terminate_macro_j++)                     \
                        {  if((x)[_terminate_macro_j] == '\n')        \
                           {  (x)[_terminate_macro_j] = '\0';         \
                              break;                                  \
                     }  }  }  while(0)
#define TERMAT(x, y) do {  int _termat_macro_j;                       \
                        for(_termat_macro_j=0;                        \
                            (x)[_termat_macro_j];                     \
                            _termat_macro_j++)                        \
                        {  if((x)[_termat_macro_j] == (y))            \
                           {  (x)[_termat_macro_j] = '\0';            \
                              break;                                  \
                     }  }  }  while(0)
#define UPPER(x) do {  int _upper_macro_i;                            \
                    for(_upper_macro_i=0;                             \
                        _upper_macro_i<(int)strlen(x) &&              \
                           (x)[_upper_macro_i];                       \
                        _upper_macro_i++)                             \
                           if(islower((x)[_upper_macro_i]))           \
                              (x)[_upper_macro_i] =                   \
                                 (char)toupper((x)[_upper_macro_i]);  \
                    }  while(0)
#define LOWER(x) do {  int _lower_macro_i;                            \
                    for(_lower_macro_i=0;                             \
                        _lower_macro_i<(int)strlen(x) &&              \
                           (x)[_lower_macro_i];                       \
                        _lower_macro_i++)                             \
                           if(isupper((x)[_lower_macro_i]))           \
                              (x)[_lower_macro_i] =                   \
                                 (char)tolower((x)[_lower_macro_i]);  \
                    }  while(0)
#define KILLLEADSPACES(y,x)                                           \
                 do \
                 {  for((y)=(x); *(y) == ' ' || *(y) == '\t'; (y)++) ; } \
                 while(0)


#define KILLTRAILSPACES(x)                                              \
do {  int _kts_macro_i;                                                 \
      _kts_macro_i = strlen(x) - 1;                                     \
      while(((x)[_kts_macro_i] == ' ' ||                                \
             (x)[_kts_macro_i] == '\t') &&                              \
            _kts_macro_i>=0)                                            \
         (_kts_macro_i)--;                                              \
      (x)[++(_kts_macro_i)] = '\0';                                     \
   }  while(0)


/* Tests for the presence of (y) in array (x) of length (l). The result
   (TRUE or FALSE) is returned in (r)
   02.08.95 Original
*/
#define TESTINARRAY(x, l, y, r)                                          \
do {                                                                     \
   int _inarray_macro_i;                                                 \
   (r) = FALSE;                                                          \
   if((x)==NULL) break;                                                  \
   for(_inarray_macro_i=0; _inarray_macro_i<(l); _inarray_macro_i++)     \
   {  if((x)[_inarray_macro_i] == (y))                                   \
      {  (r) = TRUE;                                                     \
         break;                                                          \
}  }  } while(FALSE)

/* Finds offset of item (y) in array (x) of length (l). The result
   is returned in (r) which is -1 if item not found
   02.08.95 Original
*/
#define FINDINARRAY(x, l, y, r)                                          \
do {                                                                     \
   int _inarray_macro_i;                                                 \
   (r) = (-1);                                                           \
   if((x)==NULL) break;                                                  \
   for(_inarray_macro_i=0; _inarray_macro_i<(l); _inarray_macro_i++)     \
   {  if((x)[_inarray_macro_i] == (y))                                   \
      {  (r) = _inarray_macro_i;                                         \
         break;                                                          \
}  }  } while(FALSE)


/* Used just like padterm, but doesn't touch the string if it's already
   longer than len characters
*/
#define PADMINTERM(string, len)                                \
        do {                                                   \
        if(strlen((string)) < (len)) padterm((string), (len)); \
        } while(0)

/************************************************************************/
/*>PADCHARMINTERM(string, char, length)
   ------------------------------------
   Pads a string to a specified length using char and terminates at that 
   point

   13.03.99 Original   By: ACRM
*/
#define PADCHARMINTERM(s, c, l)                                          \
do {  int _padminterm_macro_i;                                           \
      if(strlen((s)) < (l))                                              \
      {  for(_padminterm_macro_i=strlen((s));                            \
             _padminterm_macro_i<(l);                                    \
             _padminterm_macro_i++)                                      \
            (s)[_padminterm_macro_i] = (c);                              \
         (s)[(l)] = '\0';                                                \
      }  }  while(0)


/************************************************************************/
/*>DOTIFY(char *str)
   -----------------
   Macro to replace ' ' in a string with '.'

   21.04.99 Original   By: ACRM
*/
#define DOTIFY(str)                                                      \
do {                                                                     \
   char *_dotify_macro_chp;                                              \
   _dotify_macro_chp = str;                                              \
   while(*_dotify_macro_chp) {                                           \
      if(*_dotify_macro_chp==' ') *_dotify_macro_chp = '.';              \
      _dotify_macro_chp++;                                               \
}  } while(0)

/************************************************************************/
/*>DEDOTIFY(char *str)
   -------------------
   Macro to replace '.' in a string with ' '

   21.04.99 Original   By: ACRM
*/
#define DEDOTIFY(str)                                                    \
do {                                                                     \
   char *_dedotify_macro_chp;                                            \
   _dedotify_macro_chp = str;                                            \
   while(*_dedotify_macro_chp) {                                         \
      if(*_dedotify_macro_chp=='.') *_dedotify_macro_chp = ' ';          \
      _dedotify_macro_chp++;                                             \
}  } while(0)


/************************************************************************/
/*>SUBSCHAR(s, x, y)
   -----------------
   Substitute character x by character y in string s

   21.05.99 Original
*/
#define SUBSCHAR(s, x, y)                                                \
do {  char *_subschar_macro_ch = (s);                                    \
      while(*_subschar_macro_ch != '\0')                                 \
      {  if(*_subschar_macro_ch == (x)) *_subschar_macro_ch = (y);       \
         _subschar_macro_ch++;                                           \
   }  }  while(0)





/* Bit-wise operators
   02.08.95 Original
*/
#define SET(x, y)   (x) |= (y)
#define UNSET(x, y) (x) &= (~(y))
#define ISSET(x, y) ((BOOL)((x)&(y)))



#ifdef DEBUG
#define D(BUG) fprintf(stderr,"%s",BUG); fflush(stderr)
#else
#define D(BUG)
#endif


/************************** The PROMPT macro ****************************/
/* isatty() is not POSIX                                                */
#ifdef __unix
#  if defined(_POSIX_SOURCE) || !defined(_SVR4_SOURCE)
      extern int isatty(int);
#  endif
#endif

/* Default is just to print a string as a prompt                        */
#define PROMPT(in,x) printf("%s",(x))

/* More intelligent prompts for systems where we know the FILE structure*/
#ifdef __sgi
#  undef PROMPT
#  define PROMPT(in,x) do{if(isatty((in)->_file)) \
                      printf("%s",(x));}while(0)
#endif
#ifdef __linux__
#  undef PROMPT
#  define PROMPT(in,x) do{if(isatty((in)->_fileno)) \
                      printf("%s",(x));}while(0)
#endif

#endif /* _MACROS_H                                                     */

