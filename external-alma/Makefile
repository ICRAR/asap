### Notice: 
###    This Makefile assumes CASA+CASACore enviroment. 
###    For who wants to build ASAP with older CASA without CASACore 
###    environment, replace $(CASAPATH) with $(AIPSPATH). 
EXTDIR := $(shell pwd)
ATNFDIR := $(EXTDIR)/atnf

# the casa environment AIPSPATH has to be defined
CASAROOT  := $(word 1, $(CASAPATH))
CASAARCH  := $(word 2, $(CASAPATH))

COREINCD := $(shell $(shell make -s -f $(CASAROOT)/$(CASAARCH)/makedefs VARS=COREINCD eval_vars) && echo $$COREINCD)
CASALIB   := $(CASAROOT)/$(CASAARCH)/lib
USELIB := lib64

TARGETLIB   := libatnf.so
ifeq "$(CASAARCH)" "darwin"
	TARGETLIB   := libatnf.dylib
endif

# the compiler
ifndef CXX
   CXX      := g++
endif

# compiler flags
CXXFLAGS := -fPIC -O2 -g -fno-omit-frame-pointer -pipe
CXXFLAGS += -Wall -Wextra -Wcast-align -Wno-comment
# passing flags to codes for casacore
CXXOPTS := -DCASA_USECASAPATH -DCASACORE_NEEDS_RETHROW -DAIPS_STDLIB \
	-DAIPS_AUTO_STL -D_GNU_SOURCE 
# for little endian cpu (casacore + atnf flag)
CXXOPTS += -DAIPS_LITTLE_ENDIAN
# for 64-bit platforms
CXXOPTS += -DAIPS_64B
ifeq "$(CASAARCH)" "darwin"
	CXXFLAGS += -mmacosx-version-min=10.5 -dynamic -mtune=nocona \
		-mfpmath=sse -ftree-vectorize -ffast-math -fasm-blocks
	CXXOPTS += -DAIPS_DARWIN
else # so far CASAARCH is either darwin, linux_gnu, or linux_64b (Aug.2010)
	CXXFLAGS += -ansi
	CXXOPTS += -DAIPS_LINUX
# platform dependent casacore flags
	CXXOPTS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE 
endif

# casacore library dependencies
CORELIBS := -lcasa_casa -lcasa_measures -lcasa_scimath -lcasa_ms -lcasa_tables

# cfitsio
CFITSIOROOT := /usr
CFITSIOINC := -I$(CFITSIOROOT)/include/cfitsio
CFITSIOLIB := -lcfitsio

#rpfits
RPFITSROOT := /usr
RPFITSLIB := -lrpfits

# the linker
LD        := $(CXX)

LDFLAGS   := -Wl,--no-undefined -shared -Wl,-soname,$(TARGETLIB) -s -Xlinker -rpath -Xlinker $(CASALIB)
ifeq "$(CASAARCH)" "darwin"
   LDFLAGS   := -dynamiclib -Wl,-undefined -Wl,error -install_name $(CASALIB)/$(TARGETLIB) -Xlinker -rpath -Xlinker $(CASALIB)
endif

# resoleve dependencies on external package
# pks: -lcasa_casa, atnf/pks/pks_math.h
# PKSIO: -lcfitsio, -lrpfits, atnf/pks/pks_maths.h, atnf/PKSIO/****.h
#     -lcasa_casa, -lcasa_measures, -lcasa_scimath, -lcasa_ms, -lcasa_tables
INCDIRS   := -I$(COREINCD) $(CFITSIOINC) -I$(EXTDIR)
LIBS      := $(CFITSIOLIB) $(RPFITSLIB) -L$(CASALIB) $(CORELIBS)

OBJECTS := $(ATNFDIR)/pks/pks_maths.o \
	$(ATNFDIR)/PKSIO/ASTEDataset.o \
	$(ATNFDIR)/PKSIO/ASTEFXDataset.o \
	$(ATNFDIR)/PKSIO/ASTEFXReader.o \
	$(ATNFDIR)/PKSIO/ASTEReader.o \
	$(ATNFDIR)/PKSIO/FITSreader.o \
	$(ATNFDIR)/PKSIO/GBTFITSreader.o \
	$(ATNFDIR)/PKSIO/MBFITSreader.o \
	$(ATNFDIR)/PKSIO/MBrecord.o \
	$(ATNFDIR)/PKSIO/NRO45FITSReader.o \
	$(ATNFDIR)/PKSIO/NRO45Reader.o \
	$(ATNFDIR)/PKSIO/NRODataset.o \
	$(ATNFDIR)/PKSIO/NROFITSDataset.o \
	$(ATNFDIR)/PKSIO/NROOTFDataset.o \
	$(ATNFDIR)/PKSIO/NROReader.o \
	$(ATNFDIR)/PKSIO/PKSFITSreader.o \
	$(ATNFDIR)/PKSIO/PKSMS2reader.o \
	$(ATNFDIR)/PKSIO/PKSMS2writer.o \
	$(ATNFDIR)/PKSIO/PKSreader.o \
	$(ATNFDIR)/PKSIO/PKSSDwriter.o \
	$(ATNFDIR)/PKSIO/SDFITSreader.o \
	$(ATNFDIR)/PKSIO/SDFITSwriter.o #\
#	$(ATNFDIR)/PKSIO/PKSmsg.cc

HEADERS := $(ATNFDIR)/pks/pks_maths.h \
	$(ATNFDIR)/PKSIO/ASTEDataset.h \
	$(ATNFDIR)/PKSIO/ASTEFXDataset.h \
	$(ATNFDIR)/PKSIO/ASTEFXReader.h \
	$(ATNFDIR)/PKSIO/ASTEReader.h \
	$(ATNFDIR)/PKSIO/FITSreader.h \
	$(ATNFDIR)/PKSIO/GBTFITSreader.h \
	$(ATNFDIR)/PKSIO/MBFITSreader.h \
	$(ATNFDIR)/PKSIO/MBrecord.h \
	$(ATNFDIR)/PKSIO/NRO45FITSReader.h \
	$(ATNFDIR)/PKSIO/NRO45Reader.h \
	$(ATNFDIR)/PKSIO/NRODataRecord.h \
	$(ATNFDIR)/PKSIO/NRODataset.h \
	$(ATNFDIR)/PKSIO/NROFITSDataset.h \
	$(ATNFDIR)/PKSIO/NROOTFDataset.h \
	$(ATNFDIR)/PKSIO/NROReader.h \
	$(ATNFDIR)/PKSIO/PKSFITSreader.h \
	$(ATNFDIR)/PKSIO/PKSMS2reader.h \
	$(ATNFDIR)/PKSIO/PKSMS2writer.h \
	$(ATNFDIR)/PKSIO/PKSreader.h \
	$(ATNFDIR)/PKSIO/PKSrecord.h \
	$(ATNFDIR)/PKSIO/PKSSDwriter.h \
	$(ATNFDIR)/PKSIO/PKSwriter.h \
	$(ATNFDIR)/PKSIO/SDFITSreader.h \
	$(ATNFDIR)/PKSIO/SDFITSwriter.h \
	$(ATNFDIR)/PKSIO/SrcType.h #\
#	$(ATNFDIR)/PKSIO/PKSmsg.h

.SUFFIXES: .o .cc
.cc.o:
	-@ echo
	  $(CXX) $(CXXOPTS) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

.PHONY : all clean objects realclean install

all : $(TARGETLIB)

objects : $(OBJECTS)

$(TARGETLIB) : $(OBJECTS)
	-@ echo
	  $(LD) $(CXXFLAGS) $(LDFLAGS) -o $(TARGETLIB) $(OBJECTS) $(LIBS)

$(OBJECTS) : $(HEADERS)

install : 
	@if ( test -f $(TARGETLIB) ) ; then cp -f $(TARGETLIB) $(CASALIB)/. ; fi

clean :
	-rm -f $(OBJECTS)

cleanso :
	-rm -f $(TARGETLIB)

realclean : clean cleanso

