/*************************************************************************

   Program:    
   File:       Parse.h
   
   Version:    V1.8R
   Date:       11.03.94
   Function:   Include file for the command parser
   
   Copyright:  SciTech Software 1991-4
   Author:     Andrew C. R. Martin
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

************************************************************************

   Description:
   ============
   Here are defined the MAKEKEY macro, STRING and NUMBER defines, the
   KeyWd structure and return values for the parser.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  11.07.90 Original
   V1.1  08.12.92 Defines prototypes
   V1.2  16.06.93 Added memory check to MAKEKEY
   V1.3-1.6       Skipped
   V1.7  01.03.94 Added mparse()
   V1.8  11.03.94 Skipped

*************************************************************************/
#ifndef _PARSE_H
#define _PARSE_H

/************************************************************************/
/* Includes
*/
#include "MathType.h"

/************************************************************************/
/* Defines
*/

#define STRING         1    /* Defines used for the MAKEKEY macro       */
#define NUMBER         0

#define PARSE_ERRC    -1    /* Error return values from parse()         */
#define PARSE_ERRP    -2
#define PARSE_COMMENT -3

/************************************************************************/
/* Type definitions
*/
typedef struct              /* Used to store keywords for parse()       */
{
   char  *name;
   int   string, nparam;
}  KeyWd;

typedef struct              /* Used to store keywords for mparse()      */
{
   char  *name;
   int   string, minparam, maxparam;
}  MKeyWd;

/************************************************************************/
/* Macros
*/
/* Create a keyword for parse()                                         */
#define MAKEKEY(x,w,v,z) \
        (x).name = (char *)malloc((strlen(w)+2) * sizeof(char)); \
        if((x).name != NULL) strcpy((x).name,w); \
                              (x).string = v; \
                              (x).nparam = z

/* Create a keyword for mparse()                                       */
#define MAKEMKEY(x,w,v,mn,mx) \
        (x).name = (char *)malloc((strlen(w)+2) * sizeof(char)); \
        if((x).name != NULL) strcpy((x).name,w); \
                              (x).string   = v; \
                              (x).minparam = mn; \
                              (x).maxparam = mx

/************************************************************************/
/* Prototypes
*/
int parse(char *comline, int nkeys, KeyWd *keywords, REAL *REALparam,
          char **strparam);
int mparse(char *comline, int nkeys, MKeyWd *keywords, REAL *REALparam,
          char **strparam, int *nparams);
int match(char *comstring, char *string2, int *nletters);
int GetString(char *command, char *strparam);
int GetParam(char  *command, REAL *value, int *nletters);

#endif
