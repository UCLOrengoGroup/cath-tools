/*************************************************************************

   Program:    
   File:       hpgl.c
   
   Version:    V2.1R
   Date:       01.03.94
   Function:   HPGL plotting functions
   
   Copyright:  (c) SciTech Software 1991-4
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
   V1.0  25.03.91 Original
   V1.1  28.05.92 ANSIed
   V2.0  25.06.92 Modified for AMPlot2; floats->doubles
   V2.1  27.07.93 Changed some missed float->double

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "MathType.h"
#include "SysDefs.h"
#include "plotting.h"

/************************************************************************/
/* Defines and macros
*/

#define FIXVERT          /* This is used to fix the dimensions of 
                            vertical text. Not sure if this is a bug in 
                            the PLT: HPGL interpreter or is standard HPGL.
                            Used by HPGLVText().
                            Define it if it imroves your output.
                         */
#define MAXBUFF   160

/************************************************************************/
/* Globals
*/
static int     sFontHeight,               /* Font size in plotter units */
               sFontWidth;
static REAL    sFontHcm,                  /* Font size in cm            */
               sFontWcm;
static FILE    *sHPGLFile = NULL;         /* Plot file                  */

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>BOOL HPGLInit(char *filename, char *AltFont,
                 REAL xmargin,   REAL ymargin)
   --------------------------------------------
   Input:   char  *filename   HPGL file name
            char  *AltFont    Alternative font used for Greek characters
            REAL  xmargin     Unprintable x margin
            REAL  ymargin     Unprintable y margin
   Globals: REAL  PSxpicsize  X picture size
            REAL  PSypicsize  Y picture size
            REAL  PSxoffset   X offset
            REAL  PSyoffset   Y offset
   Returns: BOOL              Success

   Initialise an HPGL plot. The parameters specify the unprintable margins
   on the output device.
   25.06.92 Taken from MoG. Changed to support offsets. Added parameters.
   29.06.92 Added initialisation of alternate font.
   06.07.93 Added AltFont parameter
*/
BOOL HPGLInit(char *filename,
              char *AltFont,
              REAL xmargin, 
              REAL ymargin)
{
   char        buffer[80];
   int         xoff,
               yoff,
               HPGLAltFont;
               
   if((sHPGLFile = fopen(filename,"w")) == NULL) return(FALSE);
   
   HPGLAltFont = PS2HPGLFont(AltFont);

   xoff = (int)((PSxoffset - xmargin) * 1024);
   yoff = (int)((PSyoffset - ymargin) * 1024);

   sprintf(buffer,"IN; IP%d,%d,%d,%d;\n",xoff,yoff,
                                         (int)(PSxpicsize*1024)+xoff,
                                         (int)(PSypicsize*1024)+yoff);
   fputs(buffer,sHPGLFile);

   sprintf(buffer,"SC0,10000,0,10000;DT\\;\n");
   fputs(buffer,sHPGLFile);
   
   sprintf(buffer,"CA%d;SS;\n",HPGLAltFont);
   fputs(buffer,sHPGLFile);

   return(TRUE);
}

/************************************************************************/
/*>void HPGLPen(int num)
   ---------------------
   Input:   int  num      Pen number

   Select a Pen
   25.06.92 Taken from MoG
   08.09.92 Changed to store pen first (having seen IntroCAD output)
*/
void HPGLPen(int num)
{
   char buffer[80];
   
   sprintf(buffer,"SP;PU;SP%1d;\n",num);
   fputs(buffer,sHPGLFile);
}

/************************************************************************/
/*>void HPGLMove(REAL x, REAL y)
   -----------------------------
   Input:   REAL   x     X position (0.0--1.0)
            REAL   y     Y position (0.0--1.0)

   Move on HPGL plot
   25.06.92 Taken from MoG
*/
void HPGLMove(REAL x,
              REAL y)
{
   char buffer[80];
   
   sprintf(buffer,"PU;PA%d, %d;\n", (int)(10000*x),(int)(10000*y));
   fputs(buffer,sHPGLFile);
}

/************************************************************************/
/*>void HPGLDraw(REAL x, REAL y)
   -----------------------------
   Input:   REAL   x     X position (0.0--1.0)
            REAL   y     Y position (0.0--1.0)

   Draw on HPGL plot
   25.06.92 Taken from MoG
*/
void HPGLDraw(REAL x,
              REAL y)
{
   char buffer[80];
   
   sprintf(buffer,"PD;PA%d, %d;\n",(int)(10000*x),(int)(10000*y));
   fputs(buffer,sHPGLFile);
}

/************************************************************************/
/*>void HPGLSetDash(int style)
   ---------------------------
   Input:   int  style    Line style

   Set the line style (may be printer dependent):
         0 Solid line
         1 ............
         2 -.-.-.-.-.-.
         3 ------------
         4 -..-..-..-..
         5 --.--.--.--.

   25.06.92 Framework
   26.06.92 Original
*/
void HPGLSetDash(int style)
{
   switch(style)
   {
   case 0:                                   /* Solid line              */
      fputs("LT;\n",sHPGLFile);
      break;
   case 1:                                   /* ............            */
      fputs("LT1,2;\n",sHPGLFile);
      break;
   case 2:                                   /* -.-.-.-.-.-.            */
      fputs("LT4,3;\n",sHPGLFile);
      break;
   case 3:                                   /* ------------            */
      fputs("LT2,3;\n",sHPGLFile);
      break;
   case 4:                                   /* -..-..-..-..            */
      fputs("LT6,3;\n",sHPGLFile);
      break;
   case 5:                                   /* --.--.--.--.            */
      fputs("LT5,3;\n",sHPGLFile);
      break;
   }
}

/************************************************************************/
/*>void HPGLFont(int font, REAL size)
   ----------------------------------
   Input:   int   font       Font number
            REAL  size       Point size of font

   Set font for HPGL plot
   25.06.92 Taken from MoG
   29.06.92 Corrected CA to CS and added SS. Changed to use global width
            and height variables. Altered width to 1/2.4 * height
   27.07.93 Changed precision of floating i/o to double
*/
void HPGLFont(int    font,
              REAL   size)
{
   char   buffer[80];
   REAL width;
          
   /* Character dimensions in centimeters. We convert pts to cm, then 
      divide by 1.7 to get a better size. Width is set to half height.
   */
   sFontHcm = size * 2.54 / (1.7 * 72.0);
   sFontWcm = sFontHcm/2.4;
   
   sprintf(buffer,"PU;CS%d;SS;SI%5.3f, %5.3f;\n",
           font,sFontWcm,sFontHcm);
   fputs(buffer,sHPGLFile);

   /* Character height in scaled units                                  */
   sFontHeight = (int)((size * 10000.0) / (72.0 * PSypicsize));

   /* Character width from height in cm. This conversion is empirical!  */
   width     =  sFontHcm * 0.6154;     /* cm                            */
   width     /= 2.54;                  /* in                            */
   width     *= 10000.0/PSxpicsize;    /* Graph units                   */
   sFontWidth =  (int)width;
}

/************************************************************************/
/*>void HPGLLText(REAL x, REAL y, char *string)
   --------------------------------------------
   Input:   REAL  x        X coordinate
            REAL  y        Y coordinate
            char  *string  Text to print

   Write left justified text on HPGL plot
   25.06.92 Taken from MoG
   29.06.92 Changed to use HPGLShowText().
*/
void HPGLLText(REAL x,
               REAL y,
               char *string)
{
   char buffer[80];
   
   sprintf(buffer,"PU;PA%d, %d;",
           (int)(10000.0 * x),(int)(10000.0 * y));
   fputs(buffer,sHPGLFile);
   
   HPGLShowText(string,0,(int)(10000.0 * x),(int)(10000.0 * y));
}

/************************************************************************/
/*>void HPGLCBText(REAL x, REAL y, REAL offset, char *text)
   --------------------------------------------------------
   Input:   REAL  x        X coordinate
            REAL  y        Y coordinate
            REAL  offset   Y-offset (multiplied by font height). 
                           Move text up by this amount
            char  *string  Text to print

   Write centre-bottom justified text in HPGL

   25.06.92 Framework
   26.06.92 Original
   29.06.92 Added SimplifyText(). Changed to use HPGLShowText().
   06.07.92 Modified x-pos for wider font width
*/
void HPGLCBText(REAL x, 
                REAL y, 
                REAL offset, 
                char *text)
{
   char     buffer[80];
   int      xpos,
            ypos;
   
   xpos = (int)(10000.0 * x);
   xpos -= strlen(SimplifyText(text)) * sFontWidth / 2;
   xpos += sFontWidth / 6;
   
   ypos = (int)(10000.0 * y);
   ypos -= (int)(offset * sFontHeight);

   sprintf(buffer,"PU;PA%d, %d;", xpos, ypos);
   fputs(buffer,sHPGLFile);

   HPGLShowText(text,0,xpos,ypos);
}

/************************************************************************/
/*>void HPGLROffText(REAL x, REAL y, REAL offset, char *text)
   ----------------------------------------------------------
   Input:   REAL  x        X coordinate
            REAL  y        Y coordinate
            REAL  offset   Move left by this amount (pts)
            char  *string  Text to print

   Print right-justified text to HPGL

   25.06.92 Framework
   26.06.92 Original
   29.06.92 Added SimplifyText(). Changed to use HPGLShowText().
   06.07.92 Modified x-pos for wider font width
*/
void HPGLROffText(REAL x, 
                  REAL y, 
                  REAL offset, 
                  char *text)
{
   char     buffer[80];
   int      xpos,
            ypos;
   
   /* Base position                                                     */
   xpos = (int)(10000.0 * x);
   /* Right justify                                                     */
   xpos -= (strlen(SimplifyText(text)) * sFontWidth);
   xpos += sFontWidth/6;
   /* Convert offset from pt to plotter units                           */
   xpos += (int)((offset * 10000.0) / (72.0 * PSxpicsize));
   
   ypos = (int)(10000.0 * y);
   /* Centre y-height                                                   */
   ypos -= (int)(1.0 * sFontHeight / 3.0);
   
   sprintf(buffer,"PU;PA%d, %d;\n", xpos, ypos);
   fputs(buffer,sHPGLFile);
   HPGLShowText(text,0,xpos,ypos);
}

/************************************************************************/
/*>void HPGLLCText(REAL x, REAL y, char *text)
   -------------------------------------------
   Input:   REAL  x        X coordinate
            REAL  y        Y coordinate
            char  *text    Text to print

   Print left-centre justified text
   
   25.06.92 Framework
   26.06.92 Original
   29.06.92 Changed to use HPGLShowText().
*/
void HPGLLCText(REAL x, 
                REAL y, 
                char *text)
{
   char  buffer[80];
   int   xpos,
         ypos;
         
   xpos = (int)(10000.0 * x);
   ypos = (int)(10000.0 * y);

   /* Centre y-height                                                   */
   ypos -= (int)(1.0 * sFontHeight / 3.0);
   
   sprintf(buffer,"PU;PA%d, %d;", xpos, ypos);
   fputs(buffer,sHPGLFile);
   HPGLShowText(text,0,xpos,ypos);
}

/************************************************************************/
/*>void HPGLCTText(REAL x, REAL y, REAL offset, char *text)
   --------------------------------------------------------
   Input:   REAL  x        X coordinate
            REAL  y        Y coordinate
            REAL  offset   Y offset in points. Move text down by this.
            char  *string  Text to print

   Center Top justify text at x,y with y-offset in pts.

   25.06.92 Framework
   26.06.92 Original
   29.06.92 Added SimplifyText(). Changed to use HPGLShowText().
   01.07.92 Corrected y-positioning. Changed offset to be in pts rather
            than a multiplier of font size.
   06.07.92 Modified x-pos for wider font width
*/
void HPGLCTText(REAL x, 
                REAL y, 
                REAL offset, 
                char *text)
{
   char     buffer[80];
   int      xpos,
            ypos;
   
   xpos = (int)(10000.0 * x);
   xpos -= (strlen(SimplifyText(text)) * sFontWidth / 2.0);
   xpos += sFontWidth/6;
   
   ypos = (int)(10000.0 * y);
   /* Move down by height of font                                       */
   ypos -= sFontHeight;
   /* Move by offset                                                    */
   ypos += (int)((offset * 10000.0) / (72.0 * PSypicsize));

   sprintf(buffer,"PU;PA%d, %d;\n", xpos, ypos);
   fputs(buffer,sHPGLFile);
   
   HPGLShowText(text,0,xpos,ypos);
}

/************************************************************************/
/*>void HPGLVText(REAL x, REAL y, REAL xoff, char *text, int TitleFont,
             REAL TitleSize, char *label, int LabelFont, REAL LabelSize)
   ---------------------------------------------------------------------
   Input:   REAL x            X coordinate (in data units)
            REAL y            Y coordinate (in data units)
            REAL xoff         X-offset in pts
            char *text        Text to be written
            char *TitleFont   Font in which to write it
            REAL TitleSize    Size of font
            char *label       Label to be used to calc x offset
            char *LabelFont   Font of this label
            REAL LabelSize    Size of this label

   Write vertical text centred on x,y offset back along x by the size of
   label and by xoff in pts
   The `label' specification is used to calculate an amount by
   which to move the text back. Typically this would be the longest data
   label on the graph's Y-axis.
   The FIXVERT code is used to correct character dimensions. Not sure if
   it is a bug in the PLT: interpreter which requires this, or whether 
   it is standard HPGL. Define FIXVERT if it improves your output.

   25.06.92 Framework
   29.06.92 Original
   06.07.92 Modified x-pos for wider font width
   27.07.93 Changed precision of floating i/o to double
*/
void HPGLVText(REAL x, 
               REAL y, 
               REAL xoff, 
               char *text, 
               int  TitleFont, 
               REAL TitleSize, 
               char *label, 
               int  LabelFont, 
               REAL LabelSize)
{
   REAL   LabelWidth;
   char   buffer[240];

#ifdef FIXVERT
   REAL   width,
          fixwidth,
          height,
          fixheight;
#endif
   
   /* Find size of label                                                */
   LabelWidth  = strlen(SimplifyText(label)) * 
                 (LabelSize * 10000.0) / (2.0 * 72.0 * PSxpicsize);
   /* Convert offset from points to plotter units                       */
   xoff *= 10000.0 / (72.0 * PSxpicsize);
   
   /* Convert x & y to plotter coordinates                              */
   x *= 10000.0;
   y *= 10000.0;

   /* Modify the x-coordinate to account for the offsets                */
   x -= LabelWidth;
   x += xoff;
   
   /* Now find the y-start to centre the string vertically              */
   y -= strlen(SimplifyText(text)) * sFontWidth / 2;
   y += sFontWidth/6;
   
#ifdef FIXVERT
   /* Character dimensions in centimeters. We convert pts to cm, then 
      divide by 1.7 to get a better size. Width is set to half height
   */
   height    = TitleSize * 2.54 / (1.7 * 72.0);
   width     = height/2.0;

   fixwidth  = width  * PSxpicsize/PSypicsize;
   fixheight = height * PSypicsize/PSxpicsize;
   
   sprintf(buffer,"PU;SI%5.3f, %5.3f;\n",fixwidth,fixheight);
   fputs(buffer,sHPGLFile);
#endif

   /* Now output the text                                               */
   sprintf(buffer,"PU;PA%d,%d;DR0,1;",(int)x,(int)y);
   fputs(buffer,sHPGLFile);
   HPGLShowText(text,1,(int)x,(int)y);
   fputs("DR;\n",sHPGLFile);
   
#ifdef FIXVERT
   sprintf(buffer,"PU;SI%5.3f, %5.3f;\n",width,height);
   fputs(buffer,sHPGLFile);
#endif
}

/************************************************************************/
/*>void HPGLEnd(void)
   ------------------
   Close the HPGL plot file
   06.07.93 Original
*/
void HPGLEnd(void)
{
   fclose(sHPGLFile);
}

/************************************************************************/
/*>void HPGLShowText(char *text, BOOL orientation, int XBase, int YBase)
   ---------------------------------------------------------------------
   Input:   char  *text       Text to be displayed with control codes
            int   orientation TRUE=vertical, FALSE=horizontal
            int   XBase       Position at which to start (device coords)
            int   XBase       Position at which to start (device coords)

   Displays the text, raising or lowering as appropriate and selecting
   alternate font where required. Used by the various ...Text() routines.

   29.06.92 Original based on PostScript version.
   27.07.93 Changed precision of floating i/o to double
   11.03.94 Changed orientation to BOOL
*/
void HPGLShowText(char *text, 
                  BOOL orientation,
                  int  XBase,
                  int  YBase)
{
   char     buffer[MAXBUFF],
            OutBuff[MAXBUFF];
   int      i, j,
            chcount,
            first = TRUE;
   REAL     fixwidth,
            fixheight;
   
   /* Walk along the string                                             */
   for(i=0, j=0, chcount=0; i<strlen(text) && j<MAXBUFF-1; i++)
   {
      switch(text[i])
      {
      case '\\':           /* Should interpret next character as Greek  */
         /* Finish the current string                                   */
         if(j)
         {
            buffer[j] = '\0';
            sprintf(OutBuff,"LB%s\\;",buffer);
            fputs(OutBuff,sHPGLFile);
            j = 0;
         }
         /* Output the next character as Greek                          */
         sprintf(OutBuff,"SA;LB%c\\;SS;",text[++i]);
         fputs(OutBuff,sHPGLFile);
         chcount++;
         break;
      case '^':            /* Should raise next character               */
         /* Finish the current string                                   */
         if(j)
         {
            buffer[j] = '\0';
            sprintf(OutBuff,"LB%s\\;",buffer);
            fputs(OutBuff,sHPGLFile);
            j = 0;
         }
         
         if(first)
         {
            if(orientation) chcount++;
            first = FALSE;
         }

         /* Move to the shifted position                                */
         if(orientation)                                 /* VERTICAL    */
            sprintf(OutBuff,"PU;PA%d,%d;", XBase-sFontHeight/2,
                  (int)(YBase+(PSxpicsize*chcount*sFontWidth/PSypicsize)));
         else                                            /* HORIZONTAL  */
            sprintf(OutBuff,"PU;PA%d,%d;", XBase + chcount*sFontWidth,
                                           YBase+sFontHeight/2);
         fputs(OutBuff,sHPGLFile);
         
         /* If necessary build string                                   */
         if(text[++i] == '{')
         {
            while(text[++i] != '}' && text[i] != '\0' && j<MAXBUFF-1)
            {
               buffer[j++] = text[i];
               chcount++;
            }
         }
         else
         {
            buffer[j++] = text[i];
            chcount++;
         }
         /* Output raised string                                        */
         if(j)
         {
            buffer[j] = '\0';
            if(orientation)
            {
               fixwidth  = sFontWcm * PSxpicsize/PSypicsize;
               fixheight = sFontHcm * PSypicsize/PSxpicsize;
            }
            else
            {
               fixwidth  = sFontWcm;
               fixheight = sFontHcm;
            }
            sprintf(OutBuff,"PU;SI%f,%f;LB%s\\;PU;SI%f,%f;",
                    fixwidth,fixheight/2.0,buffer,fixwidth,fixheight);
            fputs(OutBuff,sHPGLFile);
            j = 0;
         }
         /* Reset position                                              */
         if(orientation)                                 /* VERTICAL    */
            sprintf(OutBuff,"PU;PA%d,%d;", XBase,
                  (int)(YBase+(PSxpicsize*chcount*sFontWidth/PSypicsize)));
         else                                            /* HORIZONTAL  */
            sprintf(OutBuff,"PU;PA%d,%d;", XBase + chcount*sFontWidth,
                                           YBase);
         fputs(OutBuff,sHPGLFile);
         break;
      case '_':            /* Should lower next character               */
         /* Finish the current string                                   */
         if(j)
         {
            buffer[j] = '\0';
            sprintf(OutBuff,"LB%s\\;",buffer);
            fputs(OutBuff,sHPGLFile);
            j = 0;
         }

         if(first)
         {
            if(orientation) chcount++;
            first = FALSE;
         }

         /* Move to the shifted position                                */
         if(orientation)                                 /* VERTICAL    */
            sprintf(OutBuff,"PU;PA%d,%d;", XBase+sFontHeight/4,
                  (int)(YBase+(PSxpicsize*chcount*sFontWidth/PSypicsize)));
         else                                            /* HORIZONTAL  */
            sprintf(OutBuff,"PU;PA%d,%d;", XBase + chcount*sFontWidth,
                                           YBase-sFontHeight/4);
         fputs(OutBuff,sHPGLFile);

         /* If necessary build string                                   */
         if(text[++i] == '{')
         {
            while(text[++i] != '}' && text[i] != '\0' && j<MAXBUFF-1)
            {
               buffer[j++] = text[i];
               chcount++;
            }
         }
         else
         {
            buffer[j++] = text[i];
            chcount++;
         }

         /* Output lowered string                                       */
         if(j)
         {
            buffer[j] = '\0';
            if(orientation)
            {
               fixwidth  = sFontWcm * PSxpicsize/PSypicsize;
               fixheight = sFontHcm * PSypicsize/PSxpicsize;
            }
            else
            {
               fixwidth  = sFontWcm;
               fixheight = sFontHcm;
            }
            sprintf(OutBuff,"PU;SI%f,%f;LB%s\\;PU;SI%f,%f;",
                    fixwidth,fixheight/2.0,buffer,fixwidth,fixheight);
            fputs(OutBuff,sHPGLFile);
            j = 0;
         }
         /* Reset position                                              */
         if(orientation)                                 /* VERTICAL    */
            sprintf(OutBuff,"PU;PA%d,%d;", XBase,
                 (int)(YBase+(PSxpicsize*chcount*sFontWidth/PSypicsize)));
         else                                            /* HORIZONTAL  */
            sprintf(OutBuff,"PU;PA%d,%d;", XBase + chcount*sFontWidth,
                                           YBase);
         fputs(OutBuff,sHPGLFile);

         break;
      default:             /* An ordinary character                     */
         buffer[j++] = text[i];
         chcount++;
         break;
      }
   }
   
   if(j)
   {
      buffer[j] = '\0';
      sprintf(OutBuff,"LB%s\\;",buffer);
      fputs(OutBuff,sHPGLFile);
      j = 0;
   }

   if(strlen(text)) fputs("\n",sHPGLFile);
}

