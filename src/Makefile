TARGET   := /tmp/asap/_asap.so

# the casa environment AIPSPATH has to be defined
CASAROOT  := $(word 1, $(AIPSPATH))
CASAARCH  := $(word 2, $(AIPSPATH))
CASAINC   := -I$(CASAROOT)/code/include -I$(CASAROOT)/code/casa
CASALIB   := $(CASAROOT)/$(CASAARCH)/lib


# the compiler
ifndef CXX
   CXX      := g++
endif

# compiler flags
CXXFLAGS := -fPIC -g
CXXFLAGS += -ansi -pedantic -Wno-long-long
CXXOPTS := -DAIPS_NO_TEMPLATE_SRC

# darwin specific flags
ifeq "$(CASAARCH)" "darwin"
   CXXFLAGS := -O2 -DAIPS_DARWIN
   CXXFLAGS += -ansi -pedantic -Wno-long-long -Wno-long-double
endif


# cfitsio
CFITSIOROOT := /usr/local
# If aips++ is build with LFS you need to have
# libcfitsio with lfs support
CFITSIOINC := -I$(CFITSIOROOT)/include
CFITSIOLIB := $(CFITSIOROOT)/lib/libcfitsio_lfs.a
#CFITSIOLIB := -L$(CFITSIOROOT) -lcfitsio

#rpfits
RPFITSROOT := /usr/local
RPFITSLIB := $(RPFITSROOT)/lib/librpfits.a

# This assumes all casa libs are static only (*.a)
# if not than there might by symbol resolution errors.
CASAPPLIB := $(CASALIB)/version.o \
             -L$(CASALIB) \
             -latnf -limages -lms -lcomponents -lcoordinates \
             -llattices -lfits -lmeasures -lmeasures_f \
             -ltables -lscimath -lscimath_f -lcasa  \
             $(CASALIB)/libwcs.a \
             $(RPFITSLIB) $(CFITSIOLIB) -lg2c -lstdc++

# darwin specific CASA flags
ifeq "$(CASAARCH)" "darwin"
   CASAPPLIB += -framework vecLib
else
   CASAPPLIB += -llapack -lblas
endif

# the linker
LD        := $(CXX)

LDFLAGS   := -shared -Wl$(TARGET) -s
ifeq "$(CASAARCH)" "darwin"
   LDFLAGS   := -dynamiclib -single_module
endif

# python only 2.3 has been tested
PYVERSION := 2.3
PYTHONROOT := /usr
PYTHONINC := -I$(PYTHONROOT)/include/python$(PYVERSION)
PYTHONLIB := -L$(PYTHONROOT)/lib -lpython$(PYVERSION)

# has to be build with same g++ version as casa
BOOSTROOT := /usr
BOOSTLIB  := -L$(BOOSTROOT)/lib -lboost_python
BOOSTINC  := -I$(BOOSTROOT)/include/boost


# DO NOT MODIFY AFTER THIS
INCDIRS   := -I. $(CASAINC) $(BOOSTINC) $(PYTHONINC) $(CFITSIOINC)
LIBS      := $(PYTHONLIB)  $(BOOSTLIB) $(CASAPPLIB)

CCOBJECTS := MathUtils.o \
             MathUtils2.o \
             SDLog.o \
             STAttr.o \
             SDContainer.o \
             Scantable.o \
             STFiller.o \
             STTemplates.o \
             STMath.o \
             STSubTable.o \
             STMolecules.o \
             STTcal.o \
             STFocus.o \
             STWeather.o \
             STFrequencies.o \
             STHistory.o \
             STWriter.o \
             RowAccumulator.o \
             STSelector.o \
	     SDLineFinder.o \
             SDFitter.o
PYOBJECTS := python_Scantable.o \
             python_STFiller.o \
             python_STSelector.o \
             python_STMath.o \
             python_SDFitter.o \
             python_SDLog.o \
             python_SDLineFinder.o \
	     python_asap.o

OBJECTS := $(CCOBJECTS) $(PYOBJECTS)

HEADERS   := MathUtils.h \
             SDLog.h \
             STAttr.h \
             SDContainer.h \
             Scantable.h \
             STFiller.h \
             STMath.h \
             RowAccumulator.h \
             STSubTable.h \
             STMolecules.h \
             STTcal.h \
             STFocus.h \
             STWeather.h \
             STHistory.h \
             STFrequencies.h \
             STDefs.h \
             STWriter.h \
             STSelector.h \
	     SDLineFinder.h

STATICCCLIB := libasap.a
MEMBERS := $(patsubst %,$(STATICCCLIB)(%),$(CCOBJECTS))

.SUFFIXES: .c .o .cpp .cc .cxx .C

.cpp.o:
	-@ echo
	  $(CXX) $(CXXOPTS) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cc.o:
	-@ echo
	  $(CXX) $(CXXOPTS) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cxx.o:
	-@ echo
	  $(CXX) $(CXXOPTS) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

.C.o:
	-@ echo
	  $(CXX) $(CXXOPTS) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<


.PHONY : all clean objects realclean relink

all : $(TARGET)

relink : cleanso $(TARGET)

objects : $(OBJECTS)

$(TARGET) : $(OBJECTS)
	-@ echo
	  $(LD) $(LDFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

$(OBJECTS) : $(HEADERS)

staticlib: $(MEMBERS)

clean :
	-rm -f $(OBJECTS)

cleanso :
	-rm -f $(TARGET)

realclean : clean cleanso

test: objects tasap2.o
	-@ echo
	  $(LD) -g -o /tmp/tasap2 $(OBJECTS) $(LIBS) tasap2.o
