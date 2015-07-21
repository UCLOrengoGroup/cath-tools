/*************************************************************************

   Program:    
   File:       general.h
   
   Version:    V1.12R
   Date:       30.05.02
   Function:   Header file for general purpose routines
   
   Copyright:  (c) SciTech Software 1994-2002
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      amartin@stagleys.demon.co.uk
               
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
   V1.0  11.05.94 Original    By: ACRM
   V1.1  24.08.94 Added OpenStdFiles()
   V1.2  22.09.94 Added OpenFile()
   V1.3  17.07.95 Added countchar()
   V1.4  11.09.95 Added fgetsany()
   V1.5  18.10.95 Moved YorN() to WindIO.h
   V1.6  06.11.95 Added StoreString(), InStringList() and FreeStringList()
   V1.7  15.12.95 Added QueryStrStr()
   V1.8  09.07.96 Added IndxReal()
   V1.9  18.09.96 Added padchar()
   V1.11 13.06.00 Added strcatalloc()
   V1.12 30.05.02 Added WrapString(), WrapPrint(), RightJustify(), 
                  GetWordNC() and getfield()

*************************************************************************/
#ifndef _GENERAL_H
#define _GENERAL_H

#include <stdio.h>
#include "SysDefs.h"
#include "MathType.h"

typedef struct _stringlist
{
   struct _stringlist *next;
   char               *string;
}  STRINGLIST;

void StringToLower(char *string1, char *string2);
void StringToUpper(char *string1, char *string2);
char *KillLeadSpaces(char *string);
void KillLine(FILE *fp);
void SetExtn(char *File, char *Ext);
int chindex(char *string, char ch);
void Word(char *string1, char *string2);
void WordN(char *string1, char *string2, int  MaxChar);
void padterm(char *string, int length);
void padchar(char *string, int length, char ch);
BOOL CheckExtn(char *string, char *ext);
char *ftostr(char *str, int maxlen, REAL x, int precision);

void GetFilestem(char *filename, char *stem);
int upstrcmp(char *word1, char *word2);
int upstrncmp(char *word1, char *word2, int ncomp);
char *GetWord(char *buffer, char *word, int maxsize);
BOOL OpenStdFiles(char *infile, char *outfile, FILE **in, FILE **out);
FILE *OpenFile(char *filename, char *envvar, char *mode, BOOL *noenv);
int countchar(char *string, char ch);
char *fgetsany(FILE *fp);
char *strcatalloc(char *instr, char *catstr);

STRINGLIST *StoreString(STRINGLIST *StringList, char *string);
BOOL InStringList(STRINGLIST *StringList, char *string);
void FreeStringList(STRINGLIST *StringList);

char *QueryStrStr(char *string, char *substring);

void IndexReal(REAL *arrin, int *indx, int n);

FILE *OpenOrPipe(char *filename);
int CloseOrPipe(FILE *fp);

BOOL WrapString(char *in, char *out, int maxlen);
BOOL WrapPrint(FILE *out, char *string);
void RightJustify(char *string);
char *GetWordNC(char *buffer, char *word, int maxlen);
void getfield(char *buffer, int start, int width, char *str);

#endif
