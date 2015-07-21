/*************************************************************************

   Program:    
   File:       MovePDB.c
   
   Version:    V1.2a
   Date:       06.01.11
   Function:   
   
   Copyright:  (c) SciTech Software 1993-8
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.1  01.03.94
   V1.2  27.02.98 Removed unreachable break from switch()
   V1.2a 06.01.11 Corrected description

*************************************************************************/
/* Includes
*/
#include "SysDefs.h"
#include "pdb.h"
#include "macros.h"

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
/*>BOOL MovePDB(PDB *move, PDB **from, PDB **to)
   ---------------------------------------------
   Input:   PDB    *move     PDB record to be moved
   I/O:     PDB    **from    Start of PDB linked list containing record
            PDB    **to      Start of output linked list
   Returns: BOOL             Success?

   Moves a PDB record from one linked list to another. from and to should
   point to the start of the 2 lists. If the to list hasn't been started,
   to should be NULL. Returns TRUE if moved, FALSE otherwise.

   13.05.92 Original
   19.06.92 Changed p=*to, etc. for crappy compilers
*/
BOOL MovePDB(PDB *move, PDB **from, PDB **to)
{
   PDB *p;
   BOOL ret = FALSE;
   
   if(move != NULL && *from != NULL)
   {
      /* Find the item before move in the *from list                    */
      if(move == (*from))           /* Start of list                    */
      {
         p = NULL;
      }
      else                          /* Middle of list                   */
      {
         /* Move p to item before move                                  */
         for(p = (*from); p->next && p->next != move; NEXT(p)) ;
      }
      
      /* Unlink move from *from                                         */
      if(p)          /* We're moving something in the middle of the list*/
      {
         /* Unlink move                                                 */
         p->next = move->next;
      }
      else           /* We're moving the first one in the list          */
      {
         /* If first in *from list, reset *from list                    */
         *from = move->next;
      }

      /* Add move onto the end of *to                                   */
      move->next = NULL;
      if(*to)
      {
         /* Move p to end of *to list                                   */
         for(p=(*to); p->next; NEXT(p)) ;
         /* Link in move                                                */
         p->next = move;
      }
      else
      {
         /* Initialise *to list                                         */
         *to = move;
      }
      ret = TRUE;
   }
   return(ret);
}

