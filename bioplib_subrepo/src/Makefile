# makefile for building support libraries under UNIX:

# Destination in which to place the library and the include files
LIBDEST = ${HOME}/lib
INCDEST = ${HOME}/include
SHAREDLIBDEST = ${HOME}/lib


# Define the C compiler here
CC = gcc
 
# Define any options for the C compiler here
COPT = -fPIC -ansi -Wall -pedantic -g -DGUNZIP_SUPPORT

# Define the archive/library command here
AR = ar r

# If ranlib required, set here; otherwise set to echo
RANLIB = ranlib
#RANLIB = echo

###########################################################
#          DON'T CHANGE ANYTHING BELOW THIS POINT         #
###########################################################

# Version numbers for shared libraries - strictly the MINOR numbers here
# are the minor and release numbers
# libgen
GMAJOR = 1
GMINOR = 0.1
# libbiop
BMAJOR = 1
BMINOR = 0.1

# Files for libgen.a
OFILESG = BuffInp.o CalcExtSD.o CalcSD.o CheckExtn.o CrossProd3.o \
DistPtVect.o ErrStack.o FreeStringList.o GetFilestem.o GetWord.o \
InStringList.o IndxReal.o KillLeadSpaces.o KillLine.o MatMult33_33.o \
MatMult3_33.o NComb.o NPerm.o OpenFile.o OpenStdFiles.o QueryStrStr.o \
SetExtn.o StoreString.o StringToLower.o StringToUpper.o TrueAngle.o \
VecAdd3.o VecDist.o VecLen3.o VecSub3.o Word.o WordN.o angle.o \
array2.o chindex.o countchar.o factdiv.o factorial.o fgetsany.o \
fsscanf.o ftostr.o help.o hpgl.o invert33.o openorpipe.o padchar.o \
padterm.o parse.o pearson.o pearson1.o phi.o pldist.o plotting.o \
ps.o safemem.o simpleangle.o strcatalloc.o upstrcmp.o upstrncmp.o \
WindIO.o getfield.o array3.o justify.o wrapprint.o


# Files for libbiop.a
OFILESB = AddNTerHs.o ApMatPDB.o AppendPDB.o AtomNameMatch.o \
CalcCellTrans.o CalcChiPDB.o CalcCterCoords.o CalcRMSPDB.o \
CalcTetraHCoords.o CopyPDB.o CopyPDBCoords.o CreateRotMat.o \
DNAtoAA.o DupePDB.o ExtractZonePDB.o FNam2PDB.o FindAtomInRes.o \
FindNextChainPDB.o FindNextResidue.o FindResidue.o FindResidueSpec.o \
FindZonePDB.o FitCaCbPDB.o FitCaPDB.o FitPDB.o FixCterPDB.o \
GetCGPDB.o GetCofGPDBRange.o GetCofGPDBSCRange.o GetCrystPDB.o \
GetPDBByN.o GetPDBChainLabels.o GetPDBCoor.o HAddPDB.o InPDBZone.o \
InPDBZoneSpec.o IndexPDB.o KillPDB.o KillSidechain.o KnownSeqLen.o \
LegalAtomSpec.o MovePDB.o NumericAlign.o OrderPDB.o OriginPDB.o \
PDB2Seq.o PackPDB.o ParseRes.o RdSSPDB.o RdSecPDB.o RdSeqPDB.o \
ReadCSSR.o ReadPDB.o ReadPIR.o ReadRawPIR.o ReadSimplePIR.o \
RenumAtomsPDB.o ResolPDB.o RotPDB.o SelAtPDB.o SelectCaPDB.o \
SetChi.o SetResnam.o SplitSeq.o StripHPDB.o TermPDB.o TorCoor.o \
TranslatePDB.o TrueSeqLen.o WriteCrystPDB.o WritePDB.o WrtCSSR.o \
align.o fit.o hbond.o rsc.o throne.o FitNCaCPDB.o FindRawAtomInRes.o \
WriteGromosPDB.o WholePDB.o GlyCB.o BuildAtomNeighbourPDBList.o \
FindAtomWildcardInRes.o DupeResiduePDB.o StripWatersPDB.o aalist.o \
StructurePDB.o FindHetatmResidue.o FindHetatmResidueSpec.o 


# Static libraries - the default
all : libgen.a libbiop.a 

libgen.a : $(OFILESG)
	$(AR) libgen.a $? 
	$(RANLIB) libgen.a
libbiop.a : $(OFILESB)
	$(AR) libbiop.a $? 
	$(RANLIB) libbiop.a

# Shared libraries
shared : libgens.so.$(GMAJOR).$(GMINOR) libbiops.so.$(BMAJOR).$(BMINOR)

libgens.so.$(GMAJOR).$(GMINOR) : $(OFILESG)
	$(CC) -shared -fPIC -Wl,-soname,libgens.so.$(GMAJOR) -o libgens.so.$(GMAJOR).$(GMINOR) $? -lc
libbiops.so.$(BMAJOR).$(BMINOR) : $(OFILESB)
	$(CC) -shared -fPIC -Wl,-soname,libbiops.so.$(BMAJOR) -o libbiops.so.$(BMAJOR).$(BMINOR) $? -lc

# C compilation
.c.o :
	$(CC) $(COPT) -o $@ -c $<

# Cleanup
clean :
	rm -f $(OFILESB) $(OFILESG) libgen.a libbiop.a libgens.so.$(GMAJOR).$(GMINOR) libbiops.so.$(BMAJOR).$(BMINOR)

# Installation
install :
	mkdir -p $(LIBDEST)
	mkdir -p $(INCDEST)/bioplib
	cp libgen.a $(LIBDEST)
	cp libbiop.a $(LIBDEST)
	cp *.h $(INCDEST)/bioplib

installshared : 
	mkdir -p $(SHAREDLIBDEST)
	mkdir -p $(INCDEST)/bioplib
	cp *.h $(INCDEST)/bioplib
	cp libgens.so.$(GMAJOR).$(GMINOR) $(SHAREDLIBDEST)
	cp libbiops.so.$(GMAJOR).$(GMINOR) $(SHAREDLIBDEST)
	(cd $(SHAREDLIBDEST); ln -s libgens.so.$(GMAJOR).$(GMINOR) libgens.so.$(GMAJOR); ln -s libgens.so.$(GMAJOR) libgens.so)
	(cd $(SHAREDLIBDEST); ln -s libbiops.so.$(BMAJOR).$(BMINOR) libbiops.so.$(BMAJOR); ln -s libbiops.so.$(BMAJOR) libbiops.so)

# Help
help :
	@echo "make               : Build the static libraries"
	@echo "make install       : Install the static libraries"
	@echo "make shared        : Build shared libraries"
	@echo "make installshared : Install shared libraries"
	@echo "make clean         : Remove object and library files from compilation directory"
