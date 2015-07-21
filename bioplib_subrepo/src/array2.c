/*************************************************************************

   Program:    
   File:       array.c
   
   Version:    V1.4R
   Date:       18.03.94
   Function:   Allocate and free 2D arrays
   
   Copyright:  (c) SciTech Software 1993-4
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
   Creates a 2D array where the first dimension is a set of pointers. This
   is better for passing into subroutines than the conventional C method
   of simply declaring:
      TYPE  matrix[10][10];
   which, when passed to a fuunction, loses the concept of dimensions
   unless the matrix is explicitly defined with these dimension in the
   function.
   
   This routine creates an array of pointers to 1-D arrays and can thus be
   passed to functions successfully.

**************************************************************************

   Usage:
   ======
   matrix = (TYPE **)Array2D(sizeof(TYPE), nrows, ncolumns);
   
   e.g.
   matrix = (float **)Array2D(sizeof(float), 10, 10);
   
   Returns NULL (having freed any allocated memory) if there is a problem.

**************************************************************************

   Revision History:
   =================
   V1.0  07.10.92 Original
   V1.1  29.01.93 Added includes of sysdefs.h & malloc.h for MS-DOS
   V1.2  16.06.93 Includes stdlib.h rather than malloc.h
   V1.3  01.03.94 Corrected other include file usage
   V1.4  18.03.94 Added NULL definition for systems which don't define
                  it in stdlib.h

*************************************************************************/
/* Includes
*/
#include <stdlib.h>

/************************************************************************/
/* Defines and macros
*/
#ifndef NULL
#define NULL ((void *)0)
#endif


/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>char **Array2D(int size, int dim1, int dim2)
   --------------------------------------------
   Input:   int   size    Size of an array element
            int   dim1    First dimension (number of rows)
            int   dim2    Second dimension (number of columns)
   Returns: char  **      Array of pointers. Must be cast to required 
                          type

   Create a 2D array of elements of size `size' with dimensions `dim1' 
   rows by `dim2' columns.

   07.10.92 Original
   12.07.93 Tidied and commented
*/
char **Array2D(int size, 
               int dim1, 
               int dim2)
{
   char  **array  = NULL;
   int   i;
   
   /* Allocate memory for the outer dimension array                     */
   if((array = (char **)malloc(dim1 * sizeof(char *))) == NULL)
      return(NULL);
      
   /* Set all positions to NULL                                         */
   for(i=0; i<dim1; i++)   array[i] = NULL;

   /* Allocate memory for each array in the second dimension            */
   for(i=0; i<dim1; i++)
   {
      /* If allocation fails, jump to badexit                           */
      if((array[i] = (char *)malloc(dim2 * size)) == NULL)
         goto badexit;
   }
   
   return(array);
   
badexit:
   for(i=0; i<dim1; i++)   if(array[i]) free(array[i]);
   free(array);
   return(NULL);
}

/************************************************************************/
/*>void FreeArray2D(char **array, int dim1, int dim2)
   --------------------------------------------------
   Input:   char  **    Array of pointers to be freed
            int   dim1  First dimension (number of rows)
            int   dim2  Second dimension (number of columns)

   Frees a 2D array with dimensions `dim1' rows by `dim2' columns.

   07.10.92 Original
*/
void FreeArray2D(char   **array,
                 int    dim1, 
                 int    dim2)
{
   int   i;
   
   if(array)
   {
      for(i=0; i<dim1; i++)   if(array[i]) free(array[i]);
      free(array);
   }
}
