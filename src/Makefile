TARGET   := _asap.so
ifndef CXX
   CXX      := g++
endif
CXXFLAGS := -O3 -fPIC
CXXFLAGS += -ansi
CXXFLAGS += -pedantic -Wno-long-long

CFITSIOROOT := /usr
# If aips++ is build with LFS
# CFITSIOLIB := /usr/local/lib/libcfitsio_lfs.a
CFITSIOINC := -I$(CFITSIOROOT)/include/cfitsio
CFITSIOLIB := -L$(CFITSIOROOT) -lcfitsio

CASAROOT  := $(word 1, $(AIPSPATH))
CASAINC   := -I$(CASAROOT)/code/include
CASALIB   := $(CASAROOT)/linux_gnu/lib
CASAPPLIB := $(CASALIB)/version.o \
             -L$(CASALIB) \
	     -latnf -limages -lms -lcomponents -lcoordinates \
	     -llattices -lfits -lmeasures -lmeasures_f \
	     -ltables -lscimath -lscimath_f -lcasa \
              $(CASALIB)/libwcs.a \
             -lrpfits $(CFITSIOLIB) -llapack -lblas -lg2c

LD        := $(CXX)
LDFLAGS   := -shared -Wl,--soname=$(TARGET)
LDFLAGS   += -s

BOOSTROOT := /usr
BOOSTLIB  := -L$(BOOSTROOT)/lib -lboost_python2.3
BOOSTINC  := -I$(BOOSTROOT)/include/boost

PYTHONROOT := /usr
PYTHONINC := -I$(PYTHONROOT)/include/python2.3
PYTHONLIB := -L$(PYTHONROOT)/lib -lpython2.3

INCDIRS   := -I. $(CASAINC) $(BOOSTINC) $(PYTHONINC) $(CFITSIOINC)
LIBS      := $(BOOSTLIB) $(PYTHONLIB) $(CASAPPLIB)

OBJECTS   := MathUtils.o \
             MathUtils2.o \
             SDAttr.o \
             SDContainer.o \
	     SDFitTable.o \
             SDPol.o \
             SDPol2.o \
             SDMemTable.o \
             SDReader.o \
             SDWriter.o \
             SDMath.o \
             SDFitter.o \
             SDTemplates.o \
             SDMathWrapper.o \
             SDFITSImageWriter.o \
             SDAsciiWriter.o \
	     SDLineFinder.o \
             python_SDMemTable.o \
             python_SDReader.o \
             python_SDWriter.o \
             python_SDMath.o \
             python_SDFitter.o \
             python_SDFitTable.o \
	     python_SDLineFinder.o \
             python_SD.o


HEADERS   := SDDefs.h \
             SDAttr.h \
             SDContainer.h \
	     SDFitTable.h \
             SDMemTable.h \
             SDPol.h \
             SDReader.h \
             SDWriter.h \
             SDMath.h \
             SDFitter.h \
             SDMemTableWrapper.h \
             SDMathWrapper.h \
             SDReaderWrapper.h \
             SDWriterWrapper.h \
             SDFITSImageWriter.h \
             SDAsciiWriter.h \
	     SDLineFinder.h \
             python_SD.h


.SUFFIXES: .c .o .cpp .cc .cxx .C

.cpp.o:
	-@ echo
	   $(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cc.o:
	-@ echo
	   $(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cxx.o:
	-@ echo
	   $(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

.C.o:
	-@ echo
	   $(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<


.PHONY : all clean objects realclean relink

all : $(TARGET)

relink : cleanso $(TARGET)

objects : $(OBJECTS)

$(TARGET) : $(OBJECTS)
	-@ echo
	   $(LD) $(LDFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

$(OBJECTS) : $(HEADERS)

clean :
	-rm -f $(OBJECTS)

cleanso :
	-rm -f $(TARGET)

realclean : clean cleanso
