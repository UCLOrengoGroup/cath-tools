/*************************************************************************

   Program:    
   File:       port.h
   
   Version:    V1.2
   Date:       03.04.09
   Function:   Port-specific defines to allow us to use things like 
               popen() in a clean compile
   
   Copyright:  (c) SciTech Software 1988-2009
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
   This must be included before the system includes (i.e. stdio.h etc)

**************************************************************************

   Revision History:
   =================
   V1.0  27.02.98  Original
   V1.1  17.03.09  Added Mac OS X and Windows. By: CTP
   V1.2  03.04.09  Added check for linux whether _POSIX_SOURCE already
                   defined and added further checks for MS_WINDOWS

*************************************************************************/
/***
 *** The following are necessary for getting popen() to work cleanly      
 ***/

/* Silicon graphics IRIX                                                */
#ifdef __sgi     
/*  These are for Irix6                                                 */
#   ifndef __EXTENSIONS__
#      define __EXTENSIONS__
#   endif

#   ifdef _POSIX_C_SOURCE
#      undef _POSIX_C_SOURCE
#   endif
#   define _POSIX_C_SOURCE 2

/*  These are for Irix5                                                 */
#   ifdef _POSIX_SOURCE
#      undef _POSIX_SOURCE
#   endif
#endif

/* DEC OSF/1 (Alpha)                                                    */
#ifdef __osf__
#   ifndef _XOPEN_SOURCE
#      define _XOPEN_SOURCE
#   endif
#   ifndef __STDC__
#      define __STDC__
#   endif
#endif

/* SunOS - doesn't need anything for SunOS 4.1.2                        */
#if defined(sun) && defined(sparc) && !defined(__svr4)
#endif

/* Linux                                                                */
#ifdef linux
#   ifndef _POSIX_SOURCE
#      define _POSIX_SOURCE
#   endif
#endif

/* MacOS - Doesn't need anything.                                       */
#ifdef __APPLE__
#endif

/* Windows - Does not support unix pipes.                               */
#if defined(__WIN32__) || defined(_WIN32) || defined(WIN32) ||   \
    defined(__WIN64__) || defined(_WIN64) || defined(WIN64) ||   \
    defined(__WIN95__) || defined(__NT__) || defined(__WINDOWS__) || \
    defined(msdos)     || defined(__msdos__)
#   define MS_WINDOWS 1
#   define NOPIPE
#endif
