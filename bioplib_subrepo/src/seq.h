/*************************************************************************

   Program:    
   File:       seq.h
   
   Version:    V2.10R
   Date:       27.02.07
   Function:   Header file for sequence handling
   
   Copyright:  (c) SciTech Software 1991-2007
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
   V2.0  11.03.94 Original V2 release
   V2.1  11.05.94 Added DNAtoAA() & TrueSeqLen() prototypes
   V2.2  13.05.93 Added KnownSeqLen() prototype
   V2.3  28.02.95 Added ReadRawPIR()
   V2.4  25.07.95 Added the gBioplibSeqNucleicAcid external for throne()
   V2.5  11.07.96 Added CalcMDMScore()
   V2.6  17.09.95 Added ZeroMDM()
   V2.7  26.08.97 Added macro interfaces to new DoPDB2Seq()
   V2.8  08.03.00 Added Numeric***() alignment routines
   V2.9  02.10.00 Modified DoPDB2Seq()
   V2.10 27.02.07 Added CalcMDMScoreUC() and affinealignuc()

*************************************************************************/
#ifndef _SEQ_H
#define _SEQ_H

#include "MathType.h"
#include "SysDefs.h"
#include "pdb.h"

/************************************************************************/
/* Defines and macros
 */
#define ALLOCSIZE 80  /* ReadPIR() uses this as a chunk size for 
                         allocating memory
                      */

typedef struct
{
   BOOL fragment,
        paren,
        DotInParen,
        NonExpJoin,
        UnknownPos,
        Incomplete,
        Truncated,
        Juxtapose;
   char code[16],
        name[160],
        source[160];
}  SEQINFO;

extern BOOL gBioplibSeqNucleicAcid;

#define PDB2Seq(x)         DoPDB2Seq((x), FALSE, FALSE, FALSE)
#define PDB2SeqX(x)        DoPDB2Seq((x), TRUE,  FALSE, FALSE)
#define PDB2SeqNoX(x)      DoPDB2Seq((x), FALSE, FALSE, TRUE)
#define PDB2SeqXNoX(x)     DoPDB2Seq((x), TRUE,  FALSE, TRUE)

#define PDBProt2Seq(x)     DoPDB2Seq((x), FALSE, TRUE, FALSE)
#define PDBProt2SeqX(x)    DoPDB2Seq((x), TRUE,  TRUE, FALSE)
#define PDBProt2SeqNoX(x)  DoPDB2Seq((x), FALSE, TRUE, TRUE)
#define PDBProt2SeqXNoX(x) DoPDB2Seq((x), TRUE,  TRUE, TRUE)

char throne(char *three);
char thronex(char *three);
char *onethr(char one);
char *DoPDB2Seq(PDB *pdb, BOOL DoAsxGlx, BOOL ProtOnly, BOOL NoX);
int SplitSeq(char *LinearSeq, char **seqs);
int ReadSimplePIR(FILE *fp, int  maxres, char **seqs);
int ReadPIR(FILE *fp, BOOL DoInsert, char **seqs, int maxchain, 
            SEQINFO *seqinfo, BOOL *punct, BOOL *error);
int ReadRawPIR(FILE *fp, char **seqs, int maxchain, BOOL upcase,
               SEQINFO *seqinfo, BOOL *error);
int align(char *seq1, int  length1, char *seq2, int  length2, 
          BOOL verbose, BOOL identity, int  penalty, 
          char *align1, char *align2, int  *align_len);
int affinealign(char *seq1, int  length1, char *seq2, int  length2, 
                BOOL verbose, BOOL identity, int  penalty, int penext,
                char *align1, char *align2, int  *align_len);
int CalcMDMScore(char resa, char resb);
int affinealignuc(char *seq1, int  length1, char *seq2, int  length2, 
                  BOOL verbose, BOOL identity, int  penalty, int penext,
                  char *align1, char *align2, int  *align_len);
int CalcMDMScoreUC(char resa, char resb);
BOOL ReadMDM(char *mdmfile);
int ZeroMDM(void);
char DNAtoAA(char *dna);
int TrueSeqLen(char *sequence);
int KnownSeqLen(char *sequence);
BOOL NumericReadMDM(char *mdmfile);
int NumericCalcMDMScore(int resa, int resb);
int NumericAffineAlign(int *seq1, int length1, int *seq2, int length2, 
                       BOOL verbose, BOOL identity, int penalty, int penext,
                       int *align1, int *align2, int *align_len);
#endif
