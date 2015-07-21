/*************************************************************************

   Program:    
   File:       plotting.c
   
   Version:    V1.2R
   Date:       27.02.98
   Function:   Top level HPGL/PS plotting routines
   
   Copyright:  (c) SciTech Software 1992-8
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
   These functions provide a common interface to either PostScript or
   HPGL output.
   They simplified from a set written for the Amiga which also supports
   IFF-DR2D and Amiga screen output

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  06.04.92 Original    By: ACRM
   V1.1  01.03.94 First release
   V1.2  27.02.98 Removed unreachable breaks from switch() statement

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "MathType.h"
#include "SysDefs.h"
#include "general.h"
#include "macros.h"
#include "plotting.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXTRANS     10
#define TRANS_TABLE  "HPGL.ftr"

#define MAXPEN       6
#define MAXBUFF      160

/************************************************************************/
/* Globals
*/
static struct
{
   REAL     xmin,
            ymin,
            XPScale,
            YPScale;
}  sGraph;

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>BOOL AMInitPlot(char *filename,  char *title,   int dest, 
                   REAL OutXSize,   REAL OutYSize, 
                   REAL OutXOff,    REAL OutYOff,
                   char *AltFont,   REAL xmargin,  REAL ymargin,
                   REAL DataXMin,   REAL DataYMin, 
                   REAL DataXMax,   REAL DataYMax)
   -------------------------------------------------------------
   Input:   char  *filename   File to open
            char  *title      Title for plot
            int   dest        Destination (DEST_PS or DEST_HPGL)
            REAL  OutXSize    Output plot X size (inches)
            REAL  OutYSize    Output plot Y size (inches)
            REAL  OutXOff     Output plot X offset (inches)
            REAL  OutYOff     Output plot Y offset (inches)
            char  *AltFont    Alternate font name 
            REAL  xmargin     Unprintable x margin (inches, HPGL)
                              Sensible default: 0.58
            REAL  ymargin     Unprintable y margin (inches, HPGL)
                              Sensible default: 0.1465
            REAL  DataXMin    Min data X value
            REAL  DataYMin    Min data Y value
            REAL  DataXMax    Max data X value
            REAL  DataYMax    Max data Y value
   Returns: BOOL              TRUE: OK, FALSE: Failed

   Initialise a device ready for plotting.
   
   07.05.92 Original
   25.06.92 Added HPGL support. Moved setting of PS globals to here.
   02.07.92 Put ClearWindow() in for screen plotting. Added seek to
            start of file for PS and HPGL file plotting.
   16.07.92 Added DR2D support
   17.07.92 Corrected call to InstallDR2DFonts() *before* InitDR2D().
            Changed buffer to [40]. Added bounds calc'n for DR2D.
   20.07.92 Added alternate font parameter to DR2DInit(). Added check on
            DR2DInit() return.
   22.07.92 Consider x-axis labelling in finding max x. Corrected to
            consider precision for log axes
   24.07.92 Added extras parameter to ftostr
   27.07.92 Removed the specification of EPSF offsets since the
            PSFixBoundingBox() routine takes care of all this. Increased
            size of xmax border for DR2D plots.
   06.07.93 Changed parameters
   27.02.98 Removed unreachable breaks from switch() statement
*/
BOOL AMInitPlot(char   *filename,
                char   *title,
                int    dest,
                REAL   OutXSize, 
                REAL   OutYSize, 
                REAL   OutXOff, 
                REAL   OutYOff,
                char   *AltFont,
                REAL   xmargin,
                REAL   ymargin,
                REAL   DataXMin,
                REAL   DataYMin,
                REAL   DataXMax,
                REAL   DataYMax)
{
   PSxpicsize     = OutXSize;
   PSypicsize     = OutYSize;
   PSxoffset      = OutXOff;
   PSyoffset      = OutYOff;
   
   sGraph.xmin    = DataXMin;
   sGraph.ymin    = DataYMin;
   
   sGraph.XPScale = 1.0 / (DataXMax - DataXMin);
   sGraph.YPScale = 1.0 / (DataYMax - DataYMin);
   
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      return(PSInit(filename, title, AltFont));
   case DEST_HPGL:
      return(HPGLInit(filename, AltFont, xmargin, ymargin));
   default:
      break;
   }

   return(FALSE);
}

/************************************************************************/
/*>void AMSetPen(int dest, int pen)
   --------------------------------
   Input:   int   dest      Destination
            int   pen       Pen number

   Change pen
   06.04.92 Handles screen
   07.05.92 Added PS support
   25.06.92 Added HPGL support
   16.07.92 Added DR2D support
*/   
void AMSetPen(int   dest,
              int   pen)
{
   static REAL pens[MAXPEN] = {0.5, 0.75, 1.0, 1.25, 1.5, 1.75};

   pen--;
   pen = pen % MAXPEN;

   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      PSThick(pens[pen]);
      break;
   case DEST_HPGL:
      HPGLPen(pen+1);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void AMMove(int dest, REAL x, REAL y)
   -------------------------------------
   Input:   int   dest      Destination
            REAL  x         X coordinate
            REAL  y         Y coordinate

   Move to a position specified in data coordinates.
   06.04.92 Handles screen
   10.04.92 Added log support
   29.04.92 Added check on log bounds
   07.05.92 Added PS support
   25.06.92 Added HPGL support
   16.07.92 Added DR2D support
*/
void AMMove(int  dest,
            REAL x,
            REAL y)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      PSMove(x,y);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;
      HPGLMove(x,y);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void AMDraw(int dest, REAL x, REAL y)
   -------------------------------------
   Input:   int   dest      Destination
            REAL  x         X coordinate
            REAL  y         Y coordinate

   Draw to a position specified in data coordinates.
   06.04.92 Handles screen
   10.04.92 Added log support
   29.04.92 Added check on log bounds
   07.05.92 Added PS support
   25.06.92 Added HPGL support
   16.07.92 Added DR2D support
*/
void AMDraw(int  dest,
            REAL x,
            REAL y)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;

      y = (y-sGraph.ymin) * sGraph.YPScale;

      PSDraw(x,y);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;

      y = (y-sGraph.ymin) * sGraph.YPScale;

      HPGLDraw(x,y);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void AMSetLineStyle(int dest, int style)
   ----------------------------------------
   Input:   int   dest      Destination
            int   style     Style number (0--5)

   Set the line style
   08.04.92 Framework
   07.05.92 Original (screen & PS)
   25.06.92 Added HPGL support. Removed static store of style.
   05.07.92 Corrected line patterns. They don't have commas!
   16.07.92 Added DR2D support
*/
void AMSetLineStyle(int   dest,
                    int   style)
{
   static char PSPattern[6][16] = {"",    /* PostScript line patterns   */
                                   "2",
                                   "4 1 2 1",
                                   "4",
                                   "4 3 2 2 2 3",
                                   "4 2 4 2 2 2"};
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      if(style == 0) PSClearDash();
      else           PSSetDash(PSPattern[style]);
      break;
   case DEST_HPGL:
      HPGLSetDash(style);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void AMEndLine(int dest)
   ------------------------
   Input:   int   dest      Destination

   End a line; required by PostScript actually to draw on the paper.

   07.05.92 Original
   25.06.92 Added HPGL support
   16.07.92 Added DR2D support
*/
void AMEndLine(int  dest)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      PSStroke();
      break;
   case DEST_HPGL:
      break;
   default:
      break;
   }
}


/************************************************************************/
/*>void AMSetFont(int dest, char *PSFontName, REAL FontSize)
   ---------------------------------------------------------
   Input:   int   dest         Destination
            char  *PSFontName  PostScript font name
            REAL  FontSize     Size (in points) of font

   Sets the current font using PostScript font names. If producing HPGL
   output, a lookup table is used to translate this to an HPGL font
   number

   07.04.92 Framework
   05.05.92 Original for Screen
   07.05.92 Added PS support
   25.06.92 Added HPGL support.
   29.06.92 Modified to use new PS2AmigaFont() for HPGL
   13.07.92 Added DEF_FONT parameter to SetAmigaFont()
   16.07.92 Added DR2D support
*/
void AMSetFont(int  dest, 
               char *PSFontName,
               REAL FontSize)
{
   int FontNum;
   
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      PSFont(PSFontName, FontSize);
      break;
   case DEST_HPGL:
      FontNum = PS2HPGLFont(PSFontName);
      HPGLFont(FontNum, FontSize);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void AMText(int dest, REAL x, REAL y, char *text)
   -------------------------------------------------
   Input:   int   dest      Destination
            REAL  x         X coordinate
            REAL  y         Y coordinate
            char  *text     Text to write

   Left/bottom justify text at position in data coordinates
   08.04.92 Handles screen
   10.04.92 Added log support
   29.04.92 Added check on log bounds
   07.05.92 Added PS support
   25.06.92 Added HPGL support
   16.07.92 Added DR2D support
*/
void AMText(int  dest,
            REAL x,
            REAL y,
            char *text)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
      
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      PSLText(x,y,text);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      HPGLLText(x,y,text);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void AMCBText(int dest, REAL x, REAL y, char *text)
   ---------------------------------------------------
   Input:   int   dest      Destination
            REAL  x         X coordinate
            REAL  y         Y coordinate
            char  *text     Text to print

   Centre-bottom justify text

   07.04.92 Handles screen
   10.04.92 Added log support
   29.04.92 Added check on log bounds
   07.05.92 Added PS support
   25.06.92 Added HPGL support
   16.07.92 Added DR2D support
*/
void AMCBText(int   dest,
              REAL  x,
              REAL  y,
              char  *text)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
      
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      PSCBText(x, y, 0.0, text);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      HPGLCBText(x, y, 0.0, text);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void AMRText(int dest, REAL x, REAL y, REAL offset, char *text)
   ---------------------------------------------------------------
   Input:   int   dest      Destination
            REAL  x         X coordinate
            REAL  y         Y coordinate
            REAL  offset    Move left by this amount (in points)
            char  *text     Text to print

   Right/centre justify text at position in data coordinates; offset is 
   an x-offset specified in device coordinates (pt)
   
   06.04.92 Handles screen
   07.04.92 Fix to positioning
   10.04.92 Added log support
   29.04.92 Added check on log bounds
   07.05.92 Added PS support and offset
   25.06.92 Added HPGL support
   16.07.92 Added DR2D support
*/
void AMRText(int    dest,
             REAL   x,
             REAL   y,
             REAL   offset,
             char   *text)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
      
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      PSROffText(x, y, offset, text);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      HPGLROffText(x, y, offset, text);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void AMLCText(int dest, REAL x, REAL y, char *text)
   ---------------------------------------------------
   Input:   int   dest      Destination
            REAL  x         X coordinate
            REAL  y         Y coordinate
            char  *text     Text to print

   Left/centre height justify text at position in data coordinates
   08.04.92 Handles screen
   10.04.92 Added log support
   29.04.92 Added check on log bounds
   06.05.92 Fix to height centering
   07.05.92 Added PS support
   08.05.92 Corrected Y-pos for PS
   25.06.92 Added HPGL support
   16.07.92 Added DR2D support
*/
void AMLCText(int    dest,
             REAL   x,
             REAL   y,
             char   *text)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
      
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      PSLCText(x,y,text);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      HPGLLCText(x,y,text);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void AMCTText(int dest, REAL x, REAL y, REAL CTOffset, char *text)
   ------------------------------------------------------------------
   Input:   int   dest      Destination
            REAL  x         X coordinate
            REAL  y         Y coordinate
            REAL  CTOffset  Move down by this amount (points)
            char  *text     Text to print

   Centre/top justify text at position in data coordinates. 

   06.04.92 Handles screen
   07.04.92 Fix to positioning
   10.04.92 Added log support
   29.04.92 Added check on log bounds
   07.05.92 Added PS support
   25.06.92 Added HPGL support
   01.07.92 Changed for new versions of PSCTText() and HPGLCTText() which
            take offset in points.
   16.07.92 Added DR2D support
   07.06.93 Added CTOffset param
*/
void AMCTText(int   dest,
              REAL  x,
              REAL  y,
              REAL  CTOffset,
              char  *text)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
      
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;

      y = (y-sGraph.ymin) * sGraph.YPScale;

      PSCTText(x, y, CTOffset, text);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;

      y = (y-sGraph.ymin) * sGraph.YPScale;

      HPGLCTText(x, y, CTOffset, text);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void AMEndPlot(int dest)
   ------------------------
   Input:   int   dest      Destination

   Close up a device after plotting.
   
   07.05.92 Original
   25.06.92 Added HPGL support
   01.07.92 Added blank WriteMessage() when plotting to screen
   16.07.92 Added DR2D support
*/
void AMEndPlot(int  dest)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      PSEnd();
      break;
   case DEST_HPGL:
      HPGLEnd();
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>int PS2HPGLFont(char *font)
   ---------------------------
   Input:   char  *font    PostScript font name
   Returns: int            HPGL font number

   Takes the PostScript font name and works out the best HPGL equivalent
   from a translation table. On the first call, the table is read from 
   disk and space is allocated for it. If the routine is called with a 
   NULL parameter, the space allocated for the table is freed. It is 
   quite safe to call the routine again after this has occurred; the 
   table will simply be re-read from disk.
   
   If the requested translation is unsuccessful, 0 will be returned as 
   the font number.
   
   06.07.93 Original    By: ACRM
*/
int PS2HPGLFont(char *font)
{
   int         i;
   char        buffer[MAXBUFF];
   static BOOL FirstCall   = TRUE;
   static int  NTrans      = 0;
   static struct                    /* Font translation table           */
   {
      char  *PSFont;
      int   HPGLFont;
   }  FontTable[MAXTRANS];
   
   
   /* If called with a NULL font name, we free the font table and 
      return.
   */
   if(font==NULL)
   {
      for(i=0; i<NTrans; i++)
      {
         if(FontTable[i].PSFont)
         {
            free(FontTable[i].PSFont);
            FontTable[i].PSFont = NULL;
         }
      }
      FirstCall = TRUE;
      NTrans    = 0;
      return(0);
   }
   
   
   /* If it's the first call we read the font translation table allocating
      space for the font names. Should we fail, we just drop out.
   */
   if(FirstCall)
   {
      FILE  *fp   = NULL;  /* Font translation file                     */
      
      FirstCall   = FALSE;
      
      if((fp=fopen(TRANS_TABLE,"r")) != NULL)   /* If found table       */
      {
         char  buffer[MAXBUFF];                 /* Buffer for file      */
         char  FontName[40];                    /* Font name read       */
         int   FontNum;                         /* HPGL number read     */
               
         while(fgets(buffer,MAXBUFF-1,fp))      /* Read file            */
         {
            TERMINATE(buffer);
            sscanf(buffer,"%s %d",FontName,&FontNum);
            
            if(strlen(FontName) && FontName[0] != '!')
            {
               /* We've got a font name and number. Allocate space      */
               FontTable[NTrans].PSFont = 
                  (char *)malloc((strlen(FontName)+1) * sizeof(char));
                  
               /* No room!                                              */
               if(FontTable[NTrans].PSFont == NULL) break;
               
               /* Copy in the info, down casing the font name & removing
                  leading spaces and tabs
               */
               StringToLower(FontName, buffer);
               strcpy(FontTable[NTrans].PSFont, KillLeadSpaces(buffer));
               FontTable[NTrans].HPGLFont = FontNum;
               
               /* Increment the translation count                       */
               if(++NTrans > MAXTRANS) break;
            }
         }
         fclose(fp);
      }
   }
   
   /* If we've got some translations, search for the specified font     */
   if(NTrans)
   {
      char  *ptr = NULL;
      
      StringToLower(font, buffer);
      ptr = KillLeadSpaces(buffer);

      for(i=0; i<NTrans; i++)
      {
         if(!strcmp(ptr, FontTable[i].PSFont))
         {
            return(FontTable[i].HPGLFont);
         }
      }
   }
   
   /* If no translations, or font not found, just return 0              */
   return(0);
}

/************************************************************************/
/*>char *SimplifyText(char *string)
   --------------------------------
   Input:   char  *string   String containing control codes
   Returns: char  *         String with control codes removed

   Removes control codes from a string for screen display. Also used for
   calculating string length. The returned string is stored as static
   within the routine

   06.05.92 Original
*/
char *SimplifyText(char *string)
{
   static char retstring[MAXBUFF];
   int         i, j;
   
   /* Just return the string unaltered if its too long                  */
   if(strlen(string) > MAXBUFF-1) return(string);
   
   /* Walk along the string                                             */
   for(i=0, j=0; i<strlen(string) && j<MAXBUFF-1; i++)
   {
      switch(string[i])
      {
      case '\\':           /* Should interpret next character as Greek  */
         retstring[j++] = string[++i];
         break;
      case '^':            /* Should raise next character               */
         if(string[++i] == '{')
            while(string[++i] != '}' && string[i] != '\0' && j<MAXBUFF-1)
               retstring[j++] = string[i];
         else
            retstring[j++] = string[i];
         break;
      case '_':            /* Should lower next character               */
         if(string[++i] == '{')
            while(string[++i] != '}' && string[i] != '\0' && j<MAXBUFF-1)
               retstring[j++] = string[i];
         else
            retstring[j++] = string[i];
         break;
      default:             /* An ordinary character                     */
         retstring[j++] = string[i];
         break;
      }
   }
   
   retstring[j] = '\0';

   return(retstring);
}

