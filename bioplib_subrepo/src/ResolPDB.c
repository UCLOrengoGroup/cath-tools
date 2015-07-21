/*************************************************************************

   Program:    
   File:       ResolPDB.c
   
   Version:    V1.6R
   Date:       30.05.02
   Function:   Get resolution and R-factor information out of a PDB file
   
   Copyright:  (c) UCL / Dr. Andrew C.R. Martin, 1994-2002
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
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
   See documentation for details

**************************************************************************

   Revision History:
   =================
   V1.0  28.02.94 Original
   V1.1  18.03.94 Removed extraneous printf() statement
   V1.2  17.07.96 Added check for EXPERIMENT TYPE : THEORETICAL MODEL
                  Fixed bug in searching for MODEL or NMR info
   V1.3  27.06.97 Added handing of RESOLUTION records which point you
                  to another record for the experiment type.
                  Fixed EXPERIMENT TYPE to look for NMR as well
                  Looks for EXPDTA NMR record
   V1.4  23.03.98 Added check that RESOLUTION record is in a REMARK 2
   V1.5  08.02.99 GetResolPDB() now a wrapper to GetExptl() which
                  now parses structured REMARK3 blocks and also returns
                  the Free R. Initialise some variables to 0.0
   V1.6  30.05.02 Incorporated changes from Inpharmatica - now finds
                  Electron diffraction as an experimental type. Handles
                  files without REMARK 2 correctly
   V1.7  13.12.12 Complete re-implementation of GetExptl() for remediated
                  PDB files. Old version for old PDB files available as
                  GetExptlOld()

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include "pdb.h"
#include "macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static BOOL HasText(char *ptr, char *hasWords, char *notWords);
static int SetStrucType(char *ptr);
static REAL GetNumberAfterColon(char *ptr);
static BOOL FindNextNumber(char *buffer, FILE *fp, int nlines, int nskip,
                           int ncheck, REAL *value);


/************************************************************************/
/*>BOOL GetResolPDB(FILE *fp, REAL *resolution, REAL *RFactor,
                    int *StrucType)
   -----------------------------------------------------------
   Input:   FILE *fp           PDB file pointer
   Output:  REAL *resolution   The resolution (0.0 if not applicable)
            REAL *RFactor      The R-factor (0.0 if not found)
            REAL *StrucType    Structure type:
                               STRUCTURE_TYPE_XTAL
                               STRUCTURE_TYPE_NMR
                               STRUCTURE_TYPE_MODEL
                               STRUCTURE_TYPE_UNKNOWN
   Returns: BOOL               TRUE if resolution found (even if not
                               applicable)

   This routine attempts to obtain resolution and R-factor information
   out of a PDB file. 
   It returns TRUE or FALSE to indicate whether valid information was
   found. Resolution-not-applicable structures then have the resolution
   set to zero.

   N.B.
   The resolution information returned by the routine is reliable; the
   R-factor information is stored in so many forms that it is 
   difficult to read without some form of natural language parsing, but
   we manage to handle most situations.
   The routine assumes the R-factor to be the first number after the 
   words `R-value' (or one of the other keys - see the case statement in
   the code for the valid keywords). Thus we cannot handle records of the 
   form:
   THE R-VALUE FOR 7142 REFLECTIONS BETWEEN 10.0 AND 1.97 ANGSTROMS 
   REFINEMENT CYCLE 73 IS 0.254.
   as appears in entries such as 1LZT. Here, the first number is
   the number of reflections. There is thus a kludge which sets the 
   R-factor to zero if it was read as greater than 0.5 to avoid
   this situation. In these cases, we lose the R-factor information.
   This occurs in approx 3.5% of the 1XXX PDB entries.

   25.02.94 Original   By: ACRM
   28.02.94 Added " R = " and check that VALUE wasn't B-VALUE
   17.07.96 Added check for EXPERIMENT TYPE : THEORETICAL MODEL
            Also fixed bug in searching REMARK record
   27.06.97 Added handing of RESOLUTION records which point you to 
            another record for the experiment type.
            Fixed some calls to FindNextNumber() which were checking
            an 80 character width
            Fixed EXPERIMENT TYPE to look for NMR as well
            Looks for EXPDTA NMR record
   23.03.98 Added check that RESOLUTION record is in a REMARK 2
   08.02.99 Now a wrapper to GetExptl() which also returns FreeR
*/
BOOL GetResolPDB(FILE *fp, REAL *resolution, REAL *RFactor, 
                 int *StrucType)
{
   REAL FreeR;
   
   return(GetExptl(fp, resolution, RFactor, &FreeR, StrucType));
}

/************************************************************************/
/*>BOOL GetExptl(FILE *fp, REAL *resolution, REAL *RFactor, REAL *FreeR,
                 int *StrucType)
   ---------------------------------------------------------------------
   Input:   FILE *fp           PDB file pointer
   Output:  REAL *resolution   The resolution (0.0 if not applicable)
            REAL *RFactor      The R-factor (0.0 if not found)
            REAL *FreeR        The Free R-factor (0.0 if not found)
            REAL *StrucType    Structure type:
                               STRUCTURE_TYPE_XTAL
                               STRUCTURE_TYPE_NMR
                               STRUCTURE_TYPE_MODEL
                               STRUCTURE_TYPE_UNKNOWN
   Returns: BOOL               TRUE if resolution found (even if not
                               applicable)

   This routine attempts to obtain resolution and R-factor information
   out of a PDB file. 
   It returns TRUE or FALSE to indicate whether valid information was
   found. Resolution-not-applicable structures then have the resolution
   set to zero.

   12.12.11 Original   By: ACRM
            New implementation for remediated PDB files.
            The old version is available as GetExptlOld() which handles
            old format files.
            NOTE If multiple methods specified in EXPDTA record, only the 
            first is used
            If multiple R-factors are provided in different sections, 
            then the first one is returned.
*/
BOOL GetExptl(FILE *fp, REAL *resolution, REAL *RFactor, REAL *FreeR,
              int *StrucType)
{
   char *ptr,
        buffer[MAXBUFF];

   /* Set some defaults                                                 */
   *resolution = (REAL)0.0;
   *RFactor    = (REAL)0.0;
   *FreeR      = (REAL)0.0;
   *StrucType  = STRUCTURE_TYPE_UNKNOWN;

   
   /* Make sure we're at the start of the PDB file                      */
   rewind(fp);
   
   /* Get lines from the PDB file                                       */
   while(fgets(buffer,MAXBUFF,fp))
   {
      TERMINATE(buffer);
      buffer[72] = '\0';
         
      /* Break out of the loop as soon as we hit an ATOM record         */
      if(!strncmp(buffer,"ATOM  ",6))
         break;
         
      /* See if we've found a REMARK record                             */
      if(!strncmp(buffer,"REMARK",6))
      {
         char word[80];
         int  remarkType = 0;
         
         /* See which REMARK type it is                                 */
         ptr = GetWord(buffer+6, word, 80);
         if(sscanf(word, "%d", &remarkType))
         {
            switch(remarkType)
            {
            case 2:
               if(*resolution == 0.0)
               {
                  ptr = GetWord(ptr, word, 80);
                  if(!strncmp(word, "RESOLUTION", 10))
                  {
                     ptr = GetWord(ptr, word, 80);
                     if(!sscanf(word, "%lf", resolution))
                     {
                        *resolution = 0.0;
                     }
                  }
               }
               break;
            case 3:
               if(*RFactor == 0.0)
               {
                  if(HasText(ptr, "R VALUE WORKING", "FREE"))
                  {
                     *RFactor = GetNumberAfterColon(ptr);
                  }
               }

               if(*FreeR == 0.0)
               {
                  if(HasText(ptr, "FREE R VALUE", "TEST ERROR"))
                  {
                     *FreeR = GetNumberAfterColon(ptr);
                  }
               }

               break;
            case 200:
               /* If we didn't get the structure type from EXPDTA then 
                  try here 
               */
               if(*StrucType == STRUCTURE_TYPE_UNKNOWN)
               {
                  if(HasText(ptr, "EXPERIMENT TYPE", NULL))
                  {
                     char *colon;
                     if((colon = strchr(ptr, ':'))!=NULL)
                     {
                        colon++;
                        *StrucType = SetStrucType(colon);
                     }
                  }
               }
               break;
            case 205:
               if(*StrucType == STRUCTURE_TYPE_UNKNOWN)
               {
                  *StrucType = STRUCTURE_TYPE_FIBER;
               }
               break;
            case 215:
               if(*StrucType == STRUCTURE_TYPE_UNKNOWN)
               {
                  *StrucType = STRUCTURE_TYPE_NMR;
               }
               break;
            case 217:
               if(*StrucType == STRUCTURE_TYPE_UNKNOWN)
               {
                  *StrucType = STRUCTURE_TYPE_SSNMR;
               }
               break;
            case 230:
               if(*StrucType == STRUCTURE_TYPE_UNKNOWN)
               {
                  *StrucType = STRUCTURE_TYPE_NEUTRON;
               }
               break;
            case 240:
               if(*StrucType == STRUCTURE_TYPE_UNKNOWN)
               {
                  *StrucType = STRUCTURE_TYPE_ELECTDIFF;
               }
               break;
            case 245:
               if(*StrucType == STRUCTURE_TYPE_UNKNOWN)
               {
                  *StrucType = STRUCTURE_TYPE_EM;
               }
               break;
            case 247:
               if(*StrucType == STRUCTURE_TYPE_UNKNOWN)
               {
                  *StrucType = STRUCTURE_TYPE_EM;
               }
               break;
            case 265:
               if(*StrucType == STRUCTURE_TYPE_UNKNOWN)
               {
                  *StrucType = STRUCTURE_TYPE_SOLSCAT;
               }
               break;
            }
            
         }
      }  /* End of test for it being a REMARK line                      */
      else if(!strncmp(buffer,"EXPDTA", 6))
      {
         char *semiColon;
         
         ptr = buffer+10;

         /* Terminate at semi-colon                                     */
         if((semiColon=strchr(ptr,';'))!=NULL)
         {
            *semiColon = '\0';
         }
         
         *StrucType = SetStrucType(ptr);
         
      }
      
   }  /* End of loop through PDB file                                   */

   

   /* Return successfully; the output data are already stored in the
      appropriate places
   */
   return ((*resolution > 0.0) || 
           ( *StrucType != STRUCTURE_TYPE_UNKNOWN ) );
}

/************************************************************************/
char *ReportStructureType(int StrucType)
{
   switch(StrucType)
   {
   case STRUCTURE_TYPE_UNKNOWN:
      return("Unknown");
      break;
   case STRUCTURE_TYPE_XTAL:
      return("X-ray crystal structure");
      break;
   case STRUCTURE_TYPE_NMR:
      return("NMR");
      break;
   case STRUCTURE_TYPE_MODEL:
      return("Model");
      break;
   case STRUCTURE_TYPE_ELECTDIFF:
      return("Electron Diffraction");
      break;
   case STRUCTURE_TYPE_FIBER:
      return("Fiber Diffraction");
      break;
   case STRUCTURE_TYPE_SSNMR:
      return("Solid State NMR");
      break;
   case STRUCTURE_TYPE_NEUTRON:
      return("Neutron Scattering");
      break;
   case STRUCTURE_TYPE_EM:
      return("Electron Miscroscopy");
      break;
   case STRUCTURE_TYPE_SOLSCAT:
      return("Solution Scattering");
      break;
   case STRUCTURE_TYPE_IR:
      return("Infra-red Spectroscopy");
      break;
   case STRUCTURE_TYPE_POWDER:
      return("Powder Diffraction");
      break;
   case STRUCTURE_TYPE_FRET:
      return("Fluorescence Transfer");
      break;
   default:
      return("Unknown");
      break;
   }
   return("");
}


/************************************************************************/
static REAL GetNumberAfterColon(char *ptr)
{
   char *colon;
   REAL val = 0.0;
   
   if((colon = strchr(ptr, ':'))!=NULL)
   {
      colon++;
      sscanf(colon, "%lf", &val);
   }

   return(val);
}


/************************************************************************/
static BOOL HasText(char *ptr, char *hasWords, char *notWords)
{
   char *p = ptr, 
        *h = hasWords,
        *n = notWords,
        word1[80],
        word2[80];
   int  nRequired = 0,
        nFound    = 0;
   
   
   /* Step through the words we must have                               */
   while((h=GetWord(h, word1, 80))!=NULL)
   {
      nRequired++;
      /* Step through the words in our string                           */
      p = ptr;
      while((p=GetWord(p, word2, 80))!=NULL)
      {
         if(!strcmp(word1, word2))
         {
            nFound++;
            break;
         }
      }
   }
   
   /* If we didn't find all the words return false                      */
   if(nFound != nRequired)
   {
      return(FALSE);
   }

   /* We found all the words we must have. Are 
      there any we must NOT have? 
      Step through the words we must not have 
   */
   while((n=GetWord(n, word1, 80))!=NULL)
   {
      /* Step through the words in our string                           */
      p = ptr;
      while((p=GetWord(p, word2, 80))!=NULL)
      {
         /* Return false if we have a match                             */
         if(!strcmp(word1, word2))
         {
            return(FALSE);
         }
      }
   }

   return(TRUE);
}

/************************************************************************/
static int SetStrucType(char *ptr)
{
   if(strstr(ptr, "DIFFRACTION"))
   {
      if(strstr(ptr, "X-RAY"))
      {
         return(STRUCTURE_TYPE_XTAL);
      }
      else if(strstr(ptr, "FIBER"))
      {
         return(STRUCTURE_TYPE_FIBER);
      }
      else if(strstr(ptr, "NEUTRON"))
      {
         return(STRUCTURE_TYPE_NEUTRON);
      }
      else if(strstr(ptr, "POWDER"))
      {
         return(STRUCTURE_TYPE_POWDER);
      }
   }
   else if(strstr(ptr, "ELECTRON"))
   {
      if(strstr(ptr, "CRYSTALLOGRAPHY"))
      {
         return(STRUCTURE_TYPE_ELECTDIFF);
      }
      else if(strstr(ptr, "MICROSCOPY"))
      {
         return(STRUCTURE_TYPE_EM);
      }
   }
   else if(strstr(ptr, "SOLUTION"))
   {
      if(strstr(ptr, "NMR"))
      {
         return(STRUCTURE_TYPE_NMR);
      }
      else if(strstr(ptr, "SCATTERING"))
      {
         return(STRUCTURE_TYPE_SOLSCAT);
      }
   }
   else if(strstr(ptr, "SOLID-STATE"))
   {
      if(strstr(ptr, "NMR"))
      {
         return(STRUCTURE_TYPE_SSNMR);
      }
   }
   else if(strstr(ptr, "SPECTROSCOPY"))
   {
      if(strstr(ptr, "INFRARED"))
      {
         return(STRUCTURE_TYPE_IR);
      }
   }
   else if(strstr(ptr, "FLUORESCENCE"))
   {
      if(strstr(ptr, "TRANSFER"))
      {
         return(STRUCTURE_TYPE_FRET);
      }
   }

   return(STRUCTURE_TYPE_UNKNOWN);
}


/************************************************************************/
/*>BOOL GetExptl(FILE *fp, REAL *resolution, REAL *RFactor, REAL *FreeR,
                 int *StrucType)
   ---------------------------------------------------------------------
   Input:   FILE *fp           PDB file pointer
   Output:  REAL *resolution   The resolution (0.0 if not applicable)
            REAL *RFactor      The R-factor (0.0 if not found)
            REAL *FreeR        The Free R-factor (0.0 if not found)
            REAL *StrucType    Structure type:
                               STRUCTURE_TYPE_XTAL
                               STRUCTURE_TYPE_NMR
                               STRUCTURE_TYPE_MODEL
                               STRUCTURE_TYPE_UNKNOWN
   Returns: BOOL               TRUE if resolution found (even if not
                               applicable)

   This routine attempts to obtain resolution and R-factor information
   out of a PDB file. 
   It returns TRUE or FALSE to indicate whether valid information was
   found. Resolution-not-applicable structures then have the resolution
   set to zero.

   N.B.
   The resolution information returned by the routine is reliable; the
   R-factor information is stored in so many forms that it is 
   difficult to read without some form of natural language parsing, but
   we manage to handle most situations.
   The routine assumes the R-factor to be the first number after the 
   words `R-value' (or one of the other keys - see the case statement in
   the code for the valid keywords). Thus we cannot handle records of the 
   form:
   THE R-VALUE FOR 7142 REFLECTIONS BETWEEN 10.0 AND 1.97 ANGSTROMS 
   REFINEMENT CYCLE 73 IS 0.254.
   as appears in entries such as 1LZT. Here, the first number is
   the number of reflections. There is thus a kludge which sets the 
   R-factor to zero if it was read as greater than 0.5 to avoid
   this situation. In these cases, we lose the R-factor information.
   This occurs in approx 3.5% of the 1XXX PDB entries.

   25.02.94 Original   By: ACRM
   28.02.94 Added " R = " and check that VALUE wasn't B-VALUE
   17.07.96 Added check for EXPERIMENT TYPE : THEORETICAL MODEL
            Also fixed bug in searching REMARK record
   27.06.97 Added handing of RESOLUTION records which point you to 
            another record for the experiment type.
            Fixed some calls to FindNextNumber() which were checking
            an 80 character width
            Fixed EXPERIMENT TYPE to look for NMR as well
            Looks for EXPDTA NMR record
   23.03.98 Added check that RESOLUTION record is in a REMARK 2
   08.03.99 Renamed to GetExptl() from GetResolPDB() and added
            FreeR parameter. GetResolPDB() is now a wrapper to this
            routine.
            Added additional pass which looks for the structured
            REMARK 3 records
   28.04.99 Initialise FindRefRecord et al. to zero
   18.06.99 Added other strings to the valid structured block for pass 0
            Added check for -ve R-factor
   08.09.99 Now takes the first FREE R-factor followed by 17 spaces
            rather than the last
*/
BOOL GetExptlOld(FILE *fp, REAL *resolution, REAL *RFactor, REAL *FreeR,
                 int *StrucType)
{
   BOOL ResNotApplic   = FALSE,  /* Found resolution not applicable     */
        HaveResol      = FALSE,  /* Found resolution data               */
        HaveRFac       = FALSE,  /* Found R-factor                      */
        HaveFreeR      = FALSE,  /* Found free R-factor                 */
        InAllDataBlock = FALSE;  /* Start of structured REMARK3 block   */
   
   char *ptr,
        buffer[MAXBUFF];
   int  PassNumber;

   /* Set some defaults                                                 */
   *resolution = (REAL)0.0;
   *RFactor    = (REAL)0.0;
   *FreeR      = (REAL)0.0;
   *StrucType  = STRUCTURE_TYPE_XTAL;
   
   /* We allow a series of passes of the PDB file to get the R-value 
      information. This is described in so many different ways that we 
      make checks on the most likely forms in the first pass; if that 
      fails, we make subsequent passes using less likely (and less
      definite) options.
   */
   for(PassNumber = 0; PassNumber<5; PassNumber++)
   {
      
      /* Make sure we're at the start of the PDB file                   */
      rewind(fp);
   
      /* Get lines from the PDB file                                    */
      while(fgets(buffer,MAXBUFF,fp))
      {
         TERMINATE(buffer);
         buffer[72] = '\0';
         
         /* Break out of the loop as soon as we hit an ATOM record      */
         if(!strncmp(buffer,"ATOM  ",6))
            break;
         
         /* See if we've found a REMARK record                          */
         if(!strncmp(buffer,"REMARK",6))
         {
            /* If we haven't got it already, see if this record contains
               the Resolution information
            */
            if(!HaveResol)
            {
               /* Test for a RESOLUTION sub-record not in a TITL        */
               if(((ptr = strstr(buffer,"RESOL"))!=NULL) &&
                  (strstr(buffer,"TITL ") == NULL) &&
                  (strstr(buffer,"REMARK   2")))
               {
                  HaveResol = TRUE;
                  /* If we find the word NOT, then resolution is not
                     applicable
                  */
                  if(strstr(ptr,"NOT"))
                  {
                     ResNotApplic = TRUE;
                     break;   /* Out of search through PDB              */
                  }
                  else  /* Look for the actual resolution value         */
                  {
                     if(!FindNextNumber(ptr,fp,1,10,62,resolution))
                     {
                        HaveResol = FALSE;
                     }
                  }
               }
            }
            
            /* If we've got the resolution and it is a real value,
               then start looking for the R-factor information
            */
            if(HaveResol && !ResNotApplic)
            {
               switch(PassNumber)
               {
               case 0:
                  /* 08.03.99 This pass looks for structured REMARK 3
                     records
                  */
                  if(!strncmp(buffer,"REMARK   3",10))
                  {
                     if(strstr(buffer,
                               "FIT/AGREEMENT OF MODEL WITH ALL DATA") ||
                        strstr(buffer,
                               "DATA USED IN REFINEMENT"))
                     {
                        InAllDataBlock = TRUE;
                     }
                     else if(InAllDataBlock)
                     {
                        if((ptr=strstr(buffer,
                                       "R VALUE   (WORKING + TEST"))
                           != NULL)
                        {
                           if(FindNextNumber(ptr,fp,0,50,65,RFactor))
                              HaveRFac = TRUE;
                           else
                              *RFactor = (REAL)0.0;
                        }
                        else if(!HaveRFac &&
                                ((ptr=strstr(buffer,
                                             "R VALUE          (WORKING"))
                                 != NULL))
                        {
                           if(FindNextNumber(ptr,fp,0,50,65,RFactor))
                              HaveRFac = TRUE;
                           else
                              *RFactor = (REAL)0.0;
                        }
                        
                        /* 06.09.99 Added check on HaveFreeR            */
                        if(!HaveFreeR &&
                           (ptr=strstr(buffer,
                                       "FREE R VALUE                 "))
                           != NULL)
                        {
                           if(FindNextNumber(ptr,fp,0,50,65,FreeR))
                              HaveFreeR = TRUE;
                           else
                              *FreeR = (REAL)0.0;
                        }

                        if(strstr(buffer, "NUMBER OF NON-HYDROGEN")!=NULL)
                           InAllDataBlock = FALSE;
                     }
                  }
                  break;
               case 1:
                  /* This pass we look for R and FACTOR/VALUE on the
                     same line with a space before the R
                  */
                  if((ptr=strstr(buffer," R-FAC")) != NULL || 
                     (ptr=strstr(buffer," R FAC")) != NULL || 
                     (ptr=strstr(buffer," R ="))   != NULL || 
                     (ptr=strstr(buffer," R-VAL")) != NULL || 
                     (ptr=strstr(buffer," R VAL")) != NULL)
                  {
                     if(FindNextNumber(ptr,fp,1,10,62,RFactor))
                        HaveRFac = TRUE;
                     else
                        *RFactor = (REAL)0.0;
                  }
                  break;
               case 2:
                  /* This pass we look for R and FACTOR/VALUE on the
                     same line with no space before R; beginning of
                     line, we hope.
                  */
                  if((ptr=strstr(buffer,"R-FAC")) != NULL || 
                     (ptr=strstr(buffer,"R FAC")) != NULL || 
                     (ptr=strstr(buffer,"R-VAL")) != NULL || 
                     (ptr=strstr(buffer,"R VAL")) != NULL)
                  {
                     if(FindNextNumber(ptr,fp,1,10,62,RFactor))
                        HaveRFac = TRUE;
                     else
                        *RFactor = (REAL)0.0;
                  }
                  break;
               case 3:
                  /* This pass, we look just for the word VALUE         */
                  if((ptr=strstr(buffer,"VALUE")) != NULL)
                  {
                     /* Having found VALUE, we must check that it wasn't
                        B-VALUE. Just checking the previous character
                        wasn't a - is sufficient
                     */
                     if((ptr > buffer) && (*(ptr-1) != '-'))
                     {
                        if(FindNextNumber(ptr,fp,1,10,62,RFactor))
                           HaveRFac = TRUE;
                        else
                           *RFactor = (REAL)0.0;
                     }
                  }
                  break;
               case 4:
                  /* This pass, we look just for the word FACTOR        */
                  if((ptr=strstr(buffer,"FACTOR")) != NULL)
                  {
                     if(FindNextNumber(ptr,fp,1,10,62,RFactor))
                        HaveRFac = TRUE;
                     else
                        *RFactor = (REAL)0.0;
                  }
                  break;
               }  /* End of switch                                      */

               /* If we've got an R-factor, but it's larger than 0.5,
                  then we've miss-read it. This is most likely to
                  result from the condition described in the note in
                  the header of this routine or from an error in what
                  we find in cases 2 & 3. We can't handle this
                  situation, so set value to zero.

                  18.06.99 Also check for it being negative
               */
               if(*RFactor > (REAL)0.5 || *RFactor < (REAL)0.0)
                  *RFactor = (REAL)0.0;

               /* If we've got the R-factor, break out of the search
                  through the PDB file.
               */
               if(InAllDataBlock)
               {
                  if(HaveRFac && HaveFreeR)
                     break;
               }
               else
               {
                  if(HaveRFac) 
                     break;
               }
            }  /* End of Is resolution applicable?                      */
         }  /* End of test for it being a REMARK line                   */
      }  /* End of loop through PDB file                                */

      /* After first pass, if we don't have the resolution information
         then we give up; there seems to be no interpretable information
         in the file. Output data are already set.
      */
      /* MJP 5.i.00
	 Move this to end; some XRay PDBs (e.g., 2rec) don't have a
	 remark 2 record.  This is not legal PDB, but it exists, so we
	 should parse it.
      */
      /* if(!HaveResol)
         return(FALSE); */
      
      /* Don't bother with another pass if we've found that the 
         RESOLUTION record says NOT APPLICABLE, or we've already got
         the R-factor information.
      */
      if(ResNotApplic || HaveRFac)
         break;
   }  /* End of for() each pass                                         */

   /* If it's resolution not applicable see if it's NMR or MODEL        */
   if ( ResNotApplic || !HaveResol )
   {
      REAL FindRecord    = 0.0,
           FindRefRecord = 0.0,
           ThisRecord    = 0.0;

      *StrucType = STRUCTURE_TYPE_UNKNOWN;
      
      /* Search through the same REMARK as the RESOLUTION for "MODEL",
         "NON_EXP" or "NMR" to see if it's NMR or MODEL. Also search
         all remark records for EXPERIMENT TYPE : 
      */
      /* Find the REMARK number from the buffer
         17.07.96 Corrected ptr to buffer
      */
      if(FindNextNumber(buffer, fp, 0, 0, 71, &FindRecord))
      {
         /* 27.06.97 ACRM. Also find any referenced record              */
         if(strstr(buffer+10,"REMARK"))
         {
            if(!FindNextNumber(buffer+10, fp, 0, 10, 71, &FindRefRecord))
               FindRefRecord = (-1);
         }            

         rewind(fp);
         while(fgets(buffer,MAXBUFF,fp))
         {
            if(!strncmp(buffer,"ATOM  ",6)) break;
            if(!strncmp(buffer,"EXPDTA",6))   /* 27.06.97, ACRM         */
            {
               if(strstr(buffer,"NMR"))
               {
                  *StrucType = STRUCTURE_TYPE_NMR;
                  break;
               }
	       /* 05.i.00 MJP */
	       else if (strstr (buffer,"ELECTRON MICROSCOPY"))
	       {
		   *StrucType = STRUCTURE_TYPE_ELECTDIFF;
		   break;
	       }
            }
            
            if(!strncmp(buffer,"REMARK",6))
            {
               if(FindNextNumber(buffer, fp, 0, 0, 71, &ThisRecord))
               {
                  if(((int)FindRecord    == (int)ThisRecord) ||
                     ((int)FindRefRecord == (int)ThisRecord))
                  {
                     if(strstr(buffer,"MODEL")||
                        strstr(buffer,"NON-EXP"))
                     {
                        *StrucType = STRUCTURE_TYPE_MODEL;
                        break;
                     }
                     if(strstr(buffer,"NMR"))
                     {
                        *StrucType = STRUCTURE_TYPE_NMR;
                        break;
                     }
                  }  /* In the correct REMARK block                     */
               }  /* Found a number from this REMARK record             */

               /* We also look for EXPERIMENT TYPE : MODEL / NMR
               */
               if(strstr(buffer,"EXPERIMENT"))
               {
                  if(strstr(buffer,"MODEL"))
                  {
                     *StrucType = STRUCTURE_TYPE_MODEL;
                     break;
                  }
                  else if(strstr(buffer,"NMR"))
                  {
                     *StrucType = STRUCTURE_TYPE_NMR;
                     break;
                  }
                  else if(strstr(buffer,"X-RAY"))
                  {
                     *StrucType = STRUCTURE_TYPE_XTAL;
                     break;
                  }
               }
            }  /* This is a REMARK record                               */
         }  /* while() search through PDB file                          */
      }  /* Found number in the RESOLUTION record to point to REMARK    */
   }  /* Resolution not applicable                                      */
   

   /* Return successfully; the output data are already stored in the
      appropriate places
   */
   return (HaveResol || ( *StrucType != STRUCTURE_TYPE_UNKNOWN ) );
}

/************************************************************************/
/*>static BOOL FindNextNumber(char *buffer, FILE *fp, int nlines, 
                              int nskip, int ncheck, REAL *value)
   --------------------------------------------------------------
   Find the next number which occurs in the file within nlines lines.
   First looks through the character buffer given to the routine. If this
   fails looks at the next lines in the file itself. Extra lines are
   scanned character-wise, so the whole line may not be taken out of the
   input stream. One may choose to ignore the first nskip characters from
   each new line which is read from the file. This does *not* apply to
   the initial character buffer.

   Input:   char *buffer Pointer to buffer to search first
            FILE *fp     Pointer to file to find additional lines in
            int  nlines  Number of additional lines to scan
            int  nskip   Skip this many characters at the start of
                         each new line
            int  ncheck  Check only this many characters from each new
                         line
   Output:  REAL *value  The value which we find
   Returns: BOOL         Success/failure.

   25.02.94 Original   By: ACRM
   28.02.94 Terminates at trailing decimal place
   27.06.97 Fixed potential bug when number of lines to read ahead set
            to 0 was still reading characters
   07.09.99 Added check for i > 0 in FindNextNumber when examining 
            valbuff[i-1] By: MJP
*/
static BOOL FindNextNumber(char *buffer, FILE *fp, int nlines, int nskip,
                    int ncheck, REAL *value)
{
   char *ptr,            /* Used to step through the buffers            */
        valbuff[80];     /* Buffer for copying a value into             */
   int  i,j,             /* Counter                                     */
        linecount,       /* Count lines read                            */
        chcount,         /* Count characters read from file             */
        ch;              /* Character read from file                    */
   
   linecount  = 0;
   i          = 0;
   valbuff[i] = '\0';
   
   for(ptr=buffer; *ptr; ptr++)
   {
      /* If we find a space, test the assembled string                  */
      if(*ptr == ' ' || *ptr == '\t')
      {
         valbuff[i] = '\0';
         if ( ( i > 0 ) && 
	      ( valbuff[i-1] == '.' || valbuff[i-1] == ',' ) )
            valbuff[i-1] = '\0';

         if((sscanf(valbuff,"%lf",value)) == 1)
         {
            return(TRUE);
         }
         
         
         /* If we get here, it wasn't a valid number, so reset          */
         i=0;
         continue; /* Don't bother to copy in the space                 */
      }

      /* Copy in the character                                          */
      valbuff[i++] = *ptr;
   }

   /* Our string has run out, so do the same thing, but getting 
      characters from the file instead
   */
   if(!nlines)
      return(FALSE);    /* ACRM 27.06.97                                */
   
   for(j=0; j<nskip; j++)
      fgetc(fp);
   i       = 0;
   chcount = 0;
   while(((ch = (int)fgetc(fp)) != EOF) && (linecount < nlines))
   {
      /* See if we've checked enough characters. If so, skip over all
         characters up the the next newline
      */
      if(++chcount == ncheck)
      {
         while((ch = (int)fgetc(fp)) != EOF && ch != '\n');
         if(ch == EOF) break;
      }
      
      /* If we've got a newline (including from the above skipping
         process), then reset out counters and skip over the nskip
         characters at the start of this line
      */
      if(ch == '\n')
      {
         for(j=0; j<nskip; j++)
            fgetc(fp);

         chcount = 0;
         linecount++;
      }
      

      if(ch == ' ' || ch == '\t' || ch == '\n')
      {
         valbuff[i] = '\0';
         
         if((sscanf(valbuff,"%lf",value)) == 1)
         {
            return(TRUE);
         }
         
         /* If we get here, it wasn't a valid number, so reset          */
         i=0;
         continue; /* Don't bother to copy in the space                 */
      }

      /* Copy in the character                                          */
      valbuff[i++] = ch;
   }

   /* If we get here, we failed to find a number                        */
   return(FALSE);
}


/************************************************************************/
#ifdef DEMO
int main(int argc, char **argv)
{
   FILE *fp;
   REAL resol, RFactor, FreeR;
   int  StrucType;
   
   fp = fopen(argv[1],"r");
   GetExptl(fp, &resol, &RFactor, &FreeR, &StrucType);
   
   printf("PDB:       %s\n",argv[1]);
   printf("Resol:     %f\n",resol);
   printf("RFactor:   %f\n",RFactor);
   printf("Free R:    %f\n",FreeR);
   printf("StrucType: %s\n", ReportStructureType(StrucType));
   
   return(0);
}
#endif


