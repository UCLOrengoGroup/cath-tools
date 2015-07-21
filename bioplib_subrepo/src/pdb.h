/*************************************************************************

   Program:    
   File:       pdb.h
   
   Version:    V1.47
   Date:       12.12.11
   Function:   Include file for pdb routines
   
   Copyright:  (c) SciTech Software, UCL, Reading 1993-2011
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
   V1.0  04.11.88 Original
   V1.1  22.03.90 Added Secondary structure routines
   V1.2  28.03.90 Corrected field widths for V1.2 of ReadPDB
   V1.3  04.05.90 Added clear_pdb()
   V1.4  19.06.90 Changed SEC structure to correct chain and ins widths
   V1.5  19.07.90 Added INITINDEX macro
   V1.6  09.09.91 Added define so won't screw up if included twice
   V1.7  22.09.91 Altered character sizes for alignment
   V1.8  10.06.93 Changed to use REAL rather than float. Changed 
                  order within structure
   V1.9  22.02.94 Added MAXSTDAA and MAXATINRES definitions
   V1.10 01.03.94 Added stuff for ResolPDB. Removed INIT_INDEX().
                  Added DISULPHIDE definition.
                  Added HADDINFO definition.
   V1.11 18.03.94 Added prototypes for ReadPDBOccRank() and
                  ReadPDBAtomsOccRank()
                  Added gPDBPartialOcc
   V1.12 23.05.94 Added FindNextChainPDB() prototype
   V1.13 24.08.94 Added OpenPGPFile() prototype. Added prototypes for
                  new version of FixPDB(). Added CTER styles
   V1.14 03.10.94 Added FindCofGPDBRange(), FindCofGPDBSCRange(),
                  ReadPDBALL()
   V1.15 05.10.94 Changed KillSidechain()
   V1.16 11.01.94 Added StripHPDB()
   V1.17 06.03.95 doReadPDB() is now defined here rather than static
   V1.18 17.07.95 ParseResSpec() is now a BOOL
   V1.19 24.07.95 Added FNam2PDB(), TermPDB()
   V1.20 25.07.95 Added GetPDBChainLabels()
   V1.21 08.08.95 Added FindResidueSpec() and FindNextResidue()
   V1.22 12.10.95 Added DupePDB(), CopyPDBCoords(), CalcCellTrans(),
                  GetCrystPDB(), WriteCrystPDB()
   V1.23 10.01.96 Added ExtractZonePDB()
   V1.24 08.02.96 Added FindResidue()
   V1.25 14.03.96 Added FitCaPDB(), FindAtomInRes()
   V1.26 18.06.96 Added InPDBZone() and ZONE_MODE_*. Modified prototype
                  for FindZonePDB()
   V1.27 23.07.96 Added AtomNameMatch() and LegalAtomSpec()
   V1.28 12.08.96 Added RepOneSChain() and EndRepSChain()
   V1.29 19.09.96 Added InPDBZoneSpec()
   V1.30 14.10.96 Added ReadSeqresPDB();
   V1.31 16.10.96 Added SelectCaPDB()
   V1.32 18.08.98 Changed SEC to SECSTRUC 'cos of conflict in SunOS
                  Also defines SEC macro if not defined to warn you to
                  change your code!
   V1.33 28.04.99 Added GetExptl()
   V1.34 15.02.01 Added atnam_raw[] to PDB
                  Added WriteGromosPDB(), WriteGromosPDBRecord(),
                        AtomNameRawMatch()
   V1.35 12.12.01 Added FitNCaCPDB()
   V1.36 30.05.02 Changed PDB field from 'junk' to 'record_type'
                  Added the WholePDB routines and definition
   V1.37 03.06.05 Added altpos to PDB.
                  Added altpos and atnam_raw to CLEAR_PDB
   V1.38 22.09.05 Added WritePDBRecordAtnam()
   V1.39 29.09.05 Added ParseResSpecNoUpper() and DoParseResSpec()  By: TL
   V1.40 04.01.06 Added AddCBtiGly(), AddCBtoAllGly(), 
                  StripGlyCB()      By: ACRM
   V1.41 25.01.06 Added RemoveAlternates()
   V1.42 08.11.07 Added BuildAtomNeighbourPDBList()
                        FindAtomWildcardInRes()
                        DupeResiduePDB()
   V1.43 30.04.08 Added StripWatersPDB() and ISWATER() macro
   V1.44 01.06.09 Added extras field to PDB structure
   V1.45 24.11.09 Added PDBSTRUCT, PDBCHAIN, PDBRESIDUE
                  AllocPDBStructure(), FindNextChain(),
                  FreePDBStructure()
   V1.46 26.10.11 Added FindHetatmResidueSpec() and FindHetatmResidue()
   V1.47 12.12.11 Added GetExptlOld()
                  Added ResportStructureType()
                  Added new STRUCTURE_TYPE_* defines

*************************************************************************/
#ifndef _PDB_H
#define _PDB_H

#include <stdio.h>
#include <string.h>

#include "MathType.h"
#include "SysDefs.h"
#include "general.h"

#define MAXSTDAA    21       /* Number of standard amino acids (w/ PCA)*/
#define MAXATINRES  14       /* Max number of atoms in a standard aa   */

/* This is our main PDB structure used for the PDB linked lists.

   The 'extras' field is used for flags or to attach another
   structure or array to each PDB record. For example:

   typedef struct
   {
      REAL angle;
      BOOL flag;
   }  EXTRAS;
   PDB *p;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->extras = (APTR)malloc(sizeof(EXTRAS)))==NULL)
         return(FALSE);
   }

   for(p=pdb; p!=NULL; NEXT(p))
   {
      ((EXTRAS *)p->extras)->flag = FALSE;
      ((EXTRAS *)p->extras)->angle = (REAL)0.0;
   }
*/
typedef struct pdb_entry
{
   REAL x,y,z,occ,bval;
   APTR extras;
   struct pdb_entry *next;
   int  atnum;
   int  resnum;
   char record_type[8];
   char atnam[8];
   char atnam_raw[8];
   char resnam[8];
   char insert[8];
   char chain[8];
   char altpos;
}  PDB;

typedef struct pdbresidue
{
   struct pdbresidue *next,  *prev;
   PDB               *start, *stop;
   APTR              *extras;
   int               resnum;
   char              chain[8];
   char              insert[8];
   char              resnam[8];
   char              resid[8];
} PDBRESIDUE;

typedef struct pdbchain
{
   struct pdbchain *next,  *prev;
   PDB             *start, *stop;
   PDBRESIDUE      *residues;
   APTR            *extras;
   char            chain[8];
} PDBCHAIN;

typedef struct
{
   PDB      *pdb;
   PDBCHAIN *chains;
   APTR     *extras;
} PDBSTRUCT;


#define SELECT(x,w) (x) = (char *)malloc(5 * sizeof(char)); \
                    if((x) != NULL) strncpy((x),(w),5)

typedef struct sec_entry
{
   struct sec_entry *next;
   char chain1[8];
   char ins1[8];
   char chain2[8];
   char ins2[8];
   int  res1;
   int  res2;
   char type;
}  SECSTRUC;

typedef struct _wholepdb
{
   PDB        *pdb;
   STRINGLIST *header;
   STRINGLIST *trailer;
   int        natoms;
}  WHOLEPDB;

/* This is designed to cause an error message which prints this line
   It has been tested with gcc and Irix cc and does as required in
   both cases
*/
#ifndef SEC
#   define SEC (The_type_SEC_is_now_called_SECSTRUC_You_must_change_your_code *)
#endif

typedef struct _disulphide
{
   struct _disulphide *next;
   int                res1,
                      res2;
   char               chain1[8],
                      chain2[8],
                      insert1[8],
                      insert2[8];
}  DISULPHIDE;

typedef struct
{
   int   Total,      /* Total hydrogens                                 */
         T1,         /* Type 1 C-H's                                    */
         T2,         /* Type 2 C-H2's                                   */
         T3,         /* Type 3 C-H3's                                   */
         T4,         /* Type 4 sp2 C-H's,>N-H                           */
         T5;         /* Type 5 O-H's =N-H's                             */
}  HADDINFO;

#define CLEAR_PDB(p) strcpy(p->record_type,"      "); \
                     p->atnum=0; \
                     strcpy(p->atnam,"    "); \
                     strcpy(p->atnam_raw,"    "); \
                     strcpy(p->resnam,"    "); \
                     p->resnum=0; \
                     strcpy(p->insert," "); \
                     strcpy(p->chain," "); \
                     p->x = 0.0; p->y = 0.0; p->z = 0.0; \
                     p->altpos = ' '; \
                     p->occ = 0.0; p->bval = 0.0; \
                     p->next = NULL

#define ISWATER(z)   (!strncmp((z)->resnam,"HOH",3) || \
                      !strncmp((z)->resnam,"OH2",3) || \
                      !strncmp((z)->resnam,"OHH",3) || \
                      !strncmp((z)->resnam,"DOD",3) || \
                      !strncmp((z)->resnam,"OD2",3) || \
                      !strncmp((z)->resnam,"ODD",3) || \
                      !strncmp((z)->resnam,"WAT",3))


/* These are the types returned by ResolPDB()                          */
#define STRUCTURE_TYPE_UNKNOWN   0
#define STRUCTURE_TYPE_XTAL      1
#define STRUCTURE_TYPE_NMR       2
#define STRUCTURE_TYPE_MODEL     3
#define STRUCTURE_TYPE_ELECTDIFF 4
#define STRUCTURE_TYPE_FIBER     5
#define STRUCTURE_TYPE_SSNMR     6
#define STRUCTURE_TYPE_NEUTRON   7
#define STRUCTURE_TYPE_EM        8
#define STRUCTURE_TYPE_SOLSCAT   9
#define STRUCTURE_TYPE_IR       10
#define STRUCTURE_TYPE_POWDER   11
#define STRUCTURE_TYPE_FRET     12

/* These are the styles used by FixCterPDB()                            */
#define CTER_STYLE_STD         0
#define CTER_STYLE_GROMOS      1
#define CTER_STYLE_CHARMM      2

/* Return flags from GetCrystPDB()                                      */
#define XTAL_DATA_CRYST        0x0001
#define XTAL_DATA_ORIGX        0x0002
#define XTAL_DATA_SCALE        0x0004

/* Modes for FindZonePDB()                                              */
#define ZONE_MODE_RESNUM       0
#define ZONE_MODE_SEQUENTIAL   1


/************************************************************************/
/* Globals
*/
#ifdef RSC_MAIN
   char gRSCError[80];
#else
   extern char gRSCError[80];
#endif

#ifdef READPDB_MAIN
   BOOL gPDBPartialOcc;
   BOOL gPDBMultiNMR;
#else
   extern BOOL gPDBPartialOcc;
   extern BOOL gPDBMultiNMR;
#endif

/************************************************************************/
/* Prototypes
*/
PDB *ReadPDB(FILE *fp, int *natom);
PDB *ReadPDBAll(FILE *fp, int *natom);
PDB *ReadPDBAtoms(FILE *fp, int *natom);
PDB *ReadPDBOccRank(FILE *fp, int *natom, int OccRank);
PDB *ReadPDBAtomsOccRank(FILE *fp, int *natom, int OccRank);
PDB *doReadPDB(FILE *fp, int  *natom, BOOL AllAtoms, int OccRank, 
               int ModelNum);
void WritePDB(FILE *fp, PDB *pdb);
void WritePDBRecord(FILE *fp, PDB *pdb);
void WritePDBRecordAtnam(FILE *fp, PDB  *pdb);
void WriteGromosPDB(FILE *fp, PDB *pdb);
void WriteGromosPDBRecord(FILE *fp, PDB *pdb);
void GetCofGPDB(PDB   *pdb, VEC3F *cg);
void GetCofGPDBRange(PDB *start, PDB *stop, VEC3F *cg);
void GetCofGPDBSCRange(PDB *start, PDB *stop, VEC3F *cg);
void OriginPDB(PDB *pdb);
void RotatePDB(PDB  *pdb, REAL rm[3][3]);
void TranslatePDB(PDB   *pdb, VEC3F tvect);
BOOL FitPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
BOOL FitCaPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
BOOL FitNCaCPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
BOOL FitCaCbPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
REAL CalcRMSPDB(PDB *pdb1, PDB *pdb2);
int GetPDBCoor(PDB *pdb, COOR **coor);
BOOL FindZonePDB(PDB *pdb, int start, char startinsert, int stop, char stopinsert,
                 char chain, int mode, PDB **pdb_start, PDB **pdb_stop);
int HAddPDB(FILE *fp, PDB *pdb);
int ReadPGP(FILE *fp);
FILE *OpenPGPFile(char *pgpfile, BOOL AllHyd);
PDB *SelectAtomsPDB(PDB *pdbin, int nsel, char **sel, int *natom);
PDB *StripHPDB(PDB *pdbin, int *natom);
SECSTRUC *ReadSecPDB(FILE *fp, int *nsec);
void RenumAtomsPDB(PDB *pdb);
BOOL UnPackPDB(FILE *in, FILE *out);
PDB *ReadPackedPDB(FILE *in, int *natom);
BOOL PackPDB(FILE *in, FILE *out);
void WritePackedResidue(FILE *out, PDB *start, PDB *end);
PDB *FindEndPDB(PDB *start);
PDB *FixOrderPDB(PDB *pdb, BOOL Pad, BOOL Renum);
PDB *ShuffleResPDB(PDB *start, PDB *end, BOOL Pad);
BOOL GetAtomTypes(char *resnam, char **AtomTypes);
PDB *KillPDB(PDB *pdb, PDB *prev);
void CopyPDB(PDB *out, PDB *in);
BOOL MovePDB(PDB *move, PDB **from, PDB **to);
PDB *AppendPDB(PDB *first, PDB *second);
PDB *ShuffleBB(PDB *pdb);
REAL CalcChi(PDB *pdb, int type);
PDB *GetPDBByN(PDB *pdb, int n);
void SetChi(PDB *pdb, PDB *next, REAL chi, int type);
BOOL KillSidechain(PDB *ResStart, PDB *NextRes, BOOL doCB);
void SetResnam(PDB *ResStart, PDB *NextRes, char *resnam, int resnum,   
               char *insert, char *chain);
void ApplyMatrixPDB(PDB *pdb, REAL matrix[3][3]);
BOOL GetResolPDB(FILE *fp, REAL *resolution, REAL *RFactor, 
                 int *StrucType);
BOOL GetExptl(FILE *fp, REAL *resolution, REAL *RFactor, REAL *FreeR,
              int *StrucType);
BOOL GetExptlOld(FILE *fp, REAL *resolution, REAL *RFactor, REAL *FreeR,
              int *StrucType);
char *ReportStructureType(int type);
PDB **IndexPDB(PDB *pdb, int *natom);
DISULPHIDE *ReadDisulphidesPDB(FILE *fp, BOOL *error);
BOOL ParseResSpec(char *spec, char *chain, int *resnum, char *insert);
BOOL ParseResSpecNoUpper(char *spec, char *chain, int *resnum, char *insert);
BOOL DoParseResSpec(char *spec, char *chain, int *resnum, char *insert, 
                    BOOL uppercaseresspec);
BOOL RepSChain(PDB *pdb, char *sequence, char *ChiTable, char *RefCoords);
PDB *FindNextChainPDB(PDB *pdb);
BOOL FixCterPDB(PDB *pdb, int style);
BOOL CalcCterCoords(PDB *p, PDB *ca_p, PDB *c_p, PDB *o_p);
int CalcTetraHCoords(PDB *nter, COOR *coor);
int AddNTerHs(PDB **ppdb, BOOL Charmm);
char *FNam2PDB(char *filename);
PDB *TermPDB(PDB *pdb, int length);
char *GetPDBChainLabels(PDB *pdb);
PDB *FindHetatmResidueSpec(PDB *pdb, char *resspec);
PDB *FindResidueSpec(PDB *pdb, char *resspec);
PDB *FindNextResidue(PDB *pdb);
PDB *DupePDB(PDB *in);
BOOL CopyPDBCoords(PDB *out, PDB *in);
void CalcCellTrans(VEC3F UnitCell, VEC3F CellAngles, 
                   VEC3F *xtrans, VEC3F *ytrans, VEC3F *ztrans);
int GetCrystPDB(FILE *fp, VEC3F *UnitCell, VEC3F *CellAngles,
                char *spacegroup,
                REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4]);
void WriteCrystPDB(FILE *fp, VEC3F UnitCell, VEC3F CellAngles,
                   char *spacegroup,
                   REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4]);
PDB *ExtractZonePDB(PDB *inpdb, char *chain1, int resnum1, char *insert1,
                    char *chain2, int resnum2, char *insert2);
PDB *FindResidue(PDB *pdb, char chain, int resnum, char insert);
PDB *FindHetatmResidue(PDB *pdb, char chain, int resnum, char insert);
PDB *FindAtomInRes(PDB *pdb, char *atnam);
BOOL InPDBZone(PDB *p, char chain, int resnum1, char insert1, 
               int resnum2, char insert2);
BOOL InPDBZoneSpec(PDB *p, char *resspec1, char *resspec2);
BOOL AtomNameMatch(char *atnam, char *spec, BOOL *ErrorWarn);
BOOL AtomNameRawMatch(char *atnam, char *spec, BOOL *ErrorWarn);
BOOL LegalAtomSpec(char *spec);
BOOL RepOneSChain(PDB *pdb, char *ResSpec, char aa, char *ChiTable,
                  char *RefCoords);
void EndRepSChain(void);
char **ReadSeqresPDB(FILE *fp, int *nchains);
PDB *SelectCaPDB(PDB *pdb);
char *FixAtomName(char *name, REAL occup);

void FreeWholePDB(WHOLEPDB *wpdb);
void WriteWholePDB(FILE *fp, WHOLEPDB *wpdb);
void WriteWholePDBHeader(FILE *fp, WHOLEPDB *wpdb);
void WriteWholePDBTrailer(FILE *fp, WHOLEPDB *wpdb);
WHOLEPDB *ReadWholePDB(FILE *fpin);
WHOLEPDB *ReadWholePDBAtoms(FILE *fpin);
BOOL AddCBtoGly(PDB *pdb);
BOOL AddCBtoAllGly(PDB *pdb);
PDB *StripGlyCB(PDB *pdb);
PDB *RemoveAlternates(PDB *pdb);
PDB *BuildAtomNeighbourPDBList(PDB *pdb, PDB *pRes, REAL NeighbDist);
PDB *FindAtomWildcardInRes(PDB *pdb, char *pattern);
PDB *DupeResiduePDB(PDB *in);
PDB *StripWatersPDB(PDB *pdbin, int *natom);
PDBSTRUCT *AllocPDBStructure(PDB *pdb);
PDB *FindNextChain(PDB *pdb);
void FreePDBStructure(PDBSTRUCT *pdbstruct);
#endif
