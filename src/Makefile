TARGET   := _asap.so

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
CXXFLAGS := -fPIC -O3
CXXFLAGS += -ansi -Wno-long-long -Wall
CXXOPTS := -DAIPS_NO_TEMPLATE_SRC

# darwin specific flags
ifeq "$(CASAARCH)" "darwin"
   CXXFLAGS +=-DAIPS_DARWIN
   CXXFLAGS += -Wno-long-double
endif

# cfitsio
CFITSIOROOT := /usr/local
CFITSIOINC := -I$(CFITSIOROOT)/include
CFITSIOLIB := -L$(CFITSIOROOT)/lib -lcfitsio

#rpfits
RPFITSROOT := /usr/local
RPFITSLIB := $(RPFITSROOT)/lib/librpfits.a

G2CROOT := /usr/local
G2CARCH := #$(G2CROOT)/lib/gcc/powerpc-apple-darwin8.7.0/4.2.0/libgcc.a
G2CLIB := $(G2CROOT)/lib/libg2c.a $(G2CARCH)

# This assumes all casa libs are static only (*.a)
# if not than there might by symbol resolution errors.
CASAPPLIB := $(CASALIB)/version.o \
             -L$(CASALIB) \
             -latnf -lms -lcomponents -lcoordinates \
             -llattices -lfits -lmeasures \
             -ltables -lscimath -lscimath_f -lcasa  \
             $(CASALIB)/libwcs.a \
             $(RPFITSLIB) $(CFITSIOLIB) $(G2CLIB) -lstdc++

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
   LDFLAGS   := -bundle
endif

# python
PYVERSION := 2.3
PYTHONROOT := /usr
PYTHONINC := -I$(PYTHONROOT)/include/python$(PYVERSION)
PYTHONLIB := -L$(PYTHONROOT)/lib -lpython$(PYVERSION)

# use python universal build
ifeq "$(CASAARCH)" "darwin"
   PYVERSION := 2.4
   PYTHONROOT := /Library/Frameworks/Python.framework/Versions/$(PYVERSION)
   PYTHONINC := -I$(PYTHONROOT)/include/python$(PYVERSION)
   PYTHONLIB := -L$(PYTHONROOT)/lib -framework Python
endif

# has to be build with same g++ version as casa
BOOSTROOT := /usr/local
BOOSTLIB  := -L$(BOOSTROOT)/lib -lboost_python
BOOSTINC  := -I$(BOOSTROOT)/include/boost


# DO NOT MODIFY AFTER THIS
INCDIRS   := -I. $(CASAINC)  $(CFITSIOINC) $(BOOSTINC) $(PYTHONINC)
LIBS      := $(PYTHONLIB)  $(BOOSTLIB) $(CASAPPLIB)

CCOBJECTS := MathUtils.o \
	     MathUtils2.o \
             RowAccumulator.o \
             Logger.o \
	     LineCatalog.o \
             STAttr.o \
             STHeader.o \
             STFiller.o \
             STMath.o \
             STSubTable.o \
             STMolecules.o \
             STTcal.o \
             STFocus.o \
             STWeather.o \
             STFrequencies.o \
             STHistory.o \
             STSelector.o \
	     STLineFinder.o \
             STFitter.o \
             STFitEntry.o \
             STFit.o \
             STPol.o \
             STPolLinear.o \
             STPolStokes.o \
             STWriter.o \
             STAsciiWriter.o \
             Scantable.o \
             Templates.o 
PYOBJECTS := python_Scantable.o \
             python_STFiller.o \
             python_STSelector.o \
             python_STMath.o \
             python_Fitter.o \
             python_Logger.o \
             python_STLineFinder.o \
             python_STFitEntry.o \
             python_STWriter.o \
             python_LineCatalog.o \
	     python_asap.o

OBJECTS := $(CCOBJECTS) $(PYOBJECTS)

HEADERS   := MathUtils.h \
             Logger.h \
             STAttr.h \
             STHeader.h \
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
	     STFitter.h \
	     STFit.h \
	     STFitEntry.h \
	     STLineFinder.h \
             Factory.h \
             FactoryBase.h \
	     STPol.h \
	     STPolStokes.h \
	     STPolLinear.h \
	     STWriter.h \
	     STAsciiWriter.h

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
