/*************************************************************************

   Program:    
   File:       array3.c
   
   Version:    V1.0R
   Date:       30.50.02
   Function:   Allocate and free 3D arrays
   
   Copyright:  (c) Dr. Andrew C. R. Martin, University of Reading, 2002
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
               +44 (0) 7973 431635
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************
 
   CVS Tags:
   =========
 
   Last modified by:    $Author: martin $
   Tag:                 $Name:  $
   Revision:            $Revision: 1.1 $
 
**************************************************************************

   This code is copyright solely of Dr. Andrew C. R. Martin.

**************************************************************************

   Description:
   ============
   Creates a 3D array where the first dimension is a set of pointers. This
   is better for passing into subroutines than the conventional C method
   of simply declaring:
      TYPE  matrix[10][10][10];
   which, when passed to a fuunction, loses the concept of dimensions
   unless the matrix is explicitly defined with these dimension in the
   function.
   
   This routine creates an array of pointers to 1-D arrays and can thus be
   passed to functions successfully.

**************************************************************************

   Usage:
   ======
   matrix = (TYPE ***)Array3D(sizeof(TYPE), nrows, ncolumns, nplanes);
   
   e.g.
   matrix = (float **)Array3D(sizeof(float), 10, 10, 10);
   
   Returns NULL (having freed any allocated memory) if there is a problem.

**************************************************************************

   Revision History:
   =================
   V1.0  30.05.02 Original

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
void FreeArray3D(char ***array, int dim1, int dim2, int dim3);

/************************************************************************/
/*>char ***Array3D(int size, int dim1, int dim2, int dim3)
   -------------------------------------------------------
   Input:   int   size    Size of an array element
            int   dim1    First dimension (number of rows)
            int   dim2    Second dimension (number of columns)
            int   dim3    Third dimension (number of planes)
   Returns: char  ***     Array of pointers. Must be cast to required 
                          type

   Create a 3D array of elements of size `size' with dimensions `dim1' 
   rows by `dim2' columns by `dim3' planes

   30.05.02 Original
*/
char ***Array3D(int size, int dim1, int dim2, int dim3)
{
   char  ***array  = NULL;
   int   i, j;
   
   /* Allocate memory for the outer dimension array                     */
   if((array = (char ***)malloc(dim1 * sizeof(char **))) == NULL)
      return(NULL);
      
   /* Set all positions to NULL                                         */
   for(i=0; i<dim1; i++)   array[i] = NULL;

   /* Allocate memory for each array in the second dimension            */
   for(i=0; i<dim1; i++)
   {
      /* If allocation fails, jump to badexit                           */
      if((array[i] = (char **)malloc(dim2 * sizeof(char *))) == NULL)
         goto badexit;

      /* Set all positions to NULL                                      */
      for(j=0; j<dim2; j++)   array[i][j] = NULL;

      /* Allocate memory for each array in the third dimension          */
      for(j=0; j<dim2; j++)
      {
         /* If allocation fails, jump to badexit                        */
         if((array[i][j] = (char *)malloc(dim3 * size)) == NULL)
            goto badexit;
      }
   }
   
   return(array);
   
badexit:
   FreeArray3D(array, dim1, dim2, dim3);

   return(NULL);
}

/************************************************************************/
/*>void FreeArray3D(char ***array, int dim1, int dim2, int dim3)
   -------------------------------------------------------------
   Input:   char  ***   Array of pointers to be freed
            int   dim1  First dimension (number of rows)
            int   dim2  Second dimension (number of columns)
            int   dim3  Third dimension (number of planes)

   Frees a 3D array with dimensions `dim1' rows by `dim2' columns by
   `dim3' planes.

   30.05.02 Original
*/
void FreeArray3D(char ***array, int dim1, int dim2, int dim3)
{
   int   i, j;
   
   if(array!=NULL)
   {
      for(i=0; i<dim1; i++)
      {
         if(array[i] != NULL)
         {
            for(j=0; j<dim2; j++)
            {
               if(array[i][j] != NULL)
                  free(array[i][j]);
            }
            free(array[i]);
         }
      }
      
      free(array);
   }
}
