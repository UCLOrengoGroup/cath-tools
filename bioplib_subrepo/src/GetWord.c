/*************************************************************************

   Program:    
   File:       GetWord.c
   
   Version:    V2.0
   Date:       10.06.99
   Function:   Get a space delimited word from a string
   
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
   V1.0  02.03.99 Original   By: ACRM
   V2.0  10.06.99 Complete rewrite to allow escaping of characters

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "macros.h"
#include "SysDefs.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
char *doGetWord(char *buffer, char *word, int maxlen, BOOL comma);


/************************************************************************/
/*>char *doGetWord(char *buffer, char *word, int maxlen, BOOL comma)
   -----------------------------------------------------------------
   Input:   char    *buffer     Input buffer to read words from
            int     maxlen      Max length of output word
            BOOL    comma       Treat commas like white space?
   Output:  char    *word       Word read from buffer
   Returns: char    *           Pointer to start of next word in buffer
                                or NULL

   This code is designed to be called from GetWord() or GetWordNC()

   Reads a whitespace delimted word out of buffer into word. If comma is
   TRUE, then commas are treated just like white space, otherwise they
   are treated like normal characters.

   Words containing white space may be wrapped in double inverted commas.
   A \ is used as an escape character and maybe used to escape *any*
   following character. In particular:
      "\\" -> '\'     To get a backslash
      "\ " -> ' '     To get a hard whitespace (alternatively wrap the
                      string in double inverted commas)
      "\"" -> '"'     To get a double inverted comma

   10.06.99 Original   By: ACRM (based on code from Bioplib)
*/
char *doGetWord(char *buffer, char *word, int maxlen, BOOL comma)
{
   int  i, j;
   BOOL dic    = FALSE,
        escape = FALSE;
   char *chp;
   
   /* Decrement maxlen so we can terminate correctly                    */
   maxlen--;
   
   /* Check validity of passed pointers                                 */
   if(word==NULL)
      return(NULL);

   word[0] = '\0';
   if(buffer==NULL)
      return(NULL);
   
   KILLLEADSPACES(chp, buffer);

   /* Run through each character in the input buffer                    */
   for(i=0, j=0; chp[i]; i++)
   {
      switch(chp[i])
      {
      case '\\':
         /* Use backslash as an escape character. If we've just had an
            escape, then simply store it
         */
         if(escape)
         {
            escape = FALSE;
            if(j<maxlen)
               word[j++] = chp[i];
         }
         else
         {
            escape = TRUE;
         }
         break;
      case '\"':
         /* Double inverted commas enclose strings containing white space
            If we've just had an escape then handle as a normal character,
            otherwise, toggle the dic flag
         */
         if(escape)
         {
            if(j<maxlen)
               word[j++] = chp[i];
         }
         else
         {
            TOGGLE(dic);
         }
         escape = FALSE;
         break;
      case ',':
         /* A comma is handled as white space or a normal character,
            depending on the comma flag
         */
         if(!comma)   /* Treat as default                               */
         {
            if(j<maxlen)
               word[j++] = chp[i];
            escape = FALSE;
            break;
         }
         /* Otherwise, if comma is true, just fall through to treat it
            like whitespace
         */
      case ' ':
      case '\t':
         /* If we are in double inverted commas or last char was an escape
            just handle as a normal character
         */
         if(dic || escape)
         {
            if(j<maxlen)
               word[j++] = chp[i];
         }
         else
         {
            /* Otherwise, this terminates the word, so terminate, move 
               the pointer on and return
            */
            word[j] = '\0';
            chp += i;
            KILLLEADSPACES(chp, chp);
            if(comma)
            {
               /* If we are handling commas as whitespace, then k
                  the comma if found      
               */
               if(*chp == ',') chp++;
            }
            if(*chp == '\0') chp = NULL;
            return(chp);
         }
         escape = FALSE;
         break;
      default:
         /* A normal character, copy it across                          */
         if(j<maxlen)
            word[j++] = chp[i];
         escape = FALSE;
      }
   }

   word[j] = '\0';
   return(NULL);
}

/************************************************************************/
/*>char *GetWord(char *buffer, char *word, int maxlen)
   ---------------------------------------------------
   Input:   char    *buffer     Input buffer to read words from
            int     maxlen      Max length of output word
   Output:  char    *word       Word read from buffer
   Returns: char    *           Pointer to start of next word in buffer
                                or NULL

   This code is a wrapper to doGetWord()

   Reads a whitespace/comma delimted word out of buffer into word.

   Words containing white space may be wrapped in double inverted commas.
   A \ is used as an escape character and maybe used to escape *any*
   following character. In particular:
      "\\" -> '\'     To get a backslash
      "\ " -> ' '     To get a hard whitespace (alternatively wrap the
                      string in double inverted commas)
      "\"" -> '"'     To get a double inverted comma

   10.06.99 Original   By: ACRM
*/
char *GetWord(char *buffer, char *word, int maxlen)
{
   return(doGetWord(buffer, word, maxlen, TRUE));
}

/************************************************************************/
/*>char *GetWordNC(char *buffer, char *word, int maxlen)
   -----------------------------------------------------
   Input:   char    *buffer     Input buffer to read words from
            int     maxlen      Max length of output word
   Output:  char    *word       Word read from buffer
   Returns: char    *           Pointer to start of next word in buffer
                                or NULL

   This code is a wrapper to doGetWord()

   Reads a whitespace delimted word out of buffer into word. Commas
   are treated just like normal characters.

   Words containing white space may be wrapped in double inverted commas.
   A \ is used as an escape character and maybe used to escape *any*
   following character. In particular:
      "\\" -> '\'     To get a backslash
      "\ " -> ' '     To get a hard whitespace (alternatively wrap the
                      string in double inverted commas)
      "\"" -> '"'     To get a double inverted comma

   10.06.99 Original By: ACRM
*/
char *GetWordNC(char *buffer, char *word, int maxlen)
{
   return(doGetWord(buffer, word, maxlen, FALSE));
}

