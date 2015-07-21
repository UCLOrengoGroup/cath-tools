/*************************************************************************

   Program:    
   File:       CalcExtSD.c
   
   Version:    V1.3R
   Date:       22.06.94
   Function:   Calculate mean and standard deviation
   
   Copyright:  (c) SciTech Software 1990-3
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
   This routine calculates the mean and standard deviation from a set
   of numbers. The routine is called with each value to be sampled
   and the action required is specified.

**************************************************************************

   Usage:
   ======
   #include "MathUtil.h" before using these routines

   void CalcExtSD(REAL val, int action, REAL *Sx, REAL *SxSq, 
                  int *NValues, REAL *mean, REAL *SD)
   ----------------------------------------------------------
   Input:   val     int       The value to be sampled
            action  short     0: Sample the value
                              1: Calculate & return mean and SD
                              2: Clear the sample lists
   Output:  mean    *REAL     The returned mean
            SD      *REAL     The returned standard deviation
   I/O:     Sx      *REAL     Sum of values
            SxSq    *REAL     Sum of values squared
            NValues *int      Number of values

   The output values are only set when action==1

**************************************************************************

   Revision History:
   =================
   V1.0  30.03.90 Original    By: ACRM
   V1.1  17.06.93 Modified for book
   V1.2  01.03.94 Added CalcExtSD()
   V1.3  22.06.94 Fixed for just one value

*************************************************************************/
/* Includes
*/
#include <math.h>
#include "MathType.h"

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
/*>void CalcExtSD(REAL val, int action, REAL *Sx, REAL *SxSq, 
                  int *NValues, REAL *mean, REAL *SD)
   ----------------------------------------------------------
   Calculate the mean and standard deviation from a set of numbers. 
   The routine is called with each value to be sampled and the action 
   required is specified:

   Input:   val     int       The value to be sampled
            action  short     0: Sample the value
                              1: Calculate & return mean and SD
                              2: Clear the sample lists
   Output:  mean    *REAL     The returned mean
            SD      *REAL     The returned standard deviation
   I/O:     Sx      *REAL     Sum of values
            SxSq    *REAL     Sum of values squared
            NValues *int      Number of values

   The output values are only set when action==1

   This is the same as CalcSD except that the Sx, SxSq and NValues
   variables are kept outside the function instead of being static
   within the function

   13.10.93 Original based on CalcSD   By: ACRM
   22.06.94 Fixed for only one value supplied
*/
void CalcExtSD(REAL val, int action, REAL *Sx, REAL *SxSq, 
               int *NValues, REAL *mean, REAL *SD)
{
   switch(action)
   {
   case 0:
       (*NValues)++;
       *SxSq += (val * val);
       *Sx   += val;
       break;
        
   case 1:
       *mean = *SD = (REAL)0.0;
       if(*NValues > 0)
          *mean = (*Sx) / (*NValues);
       if(*NValues > 1)
          *SD   = sqrt((*SxSq - ((*Sx) * (*Sx)) / (*NValues)) /
                       (*NValues - 1));
       break;
        
   case 2:
       *SxSq = 0.0;
       *Sx   = 0.0;
       *NValues = 0;
   }
}

