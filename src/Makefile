TARGET   := /tmp/_asap.so
ifndef CXX
   CXX      := g++
endif
CXXFLAGS := -O3 -fPIC
CXXFLAGS += -ansi
CXXFLAGS += -pedantic -Wno-long-long
CXXOPTS := -DAIPS_NO_TEMPLATE_SRC

CFITSIOROOT := /usr/local If aips++ is build with LFS you need to have
# libcfitsio with lfs support
CFITSIOLIB := /usr/local/lib/libcfitsio_lfs.a
#CFITSIOINC := -I$(CFITSIOROOT)/include/cfitsio
#CFITSIOLIB := -L$(CFITSIOROOT) -lcfitsio


CASAROOT  := $(word 1, $(AIPSPATH))
CASAARCH  := $(word 2, $(AIPSPATH))
CASAINC   := -I$(CASAROOT)/code/include -I$(CASAROOT)/code/casa/wcslib
CASALIB   := $(CASAROOT)/$(CASAARCH)/lib


# This assumes all casa libs are static only (*.a)
# if not than there might by symbol resolution errors.
CASAPPLIB := $(CASALIB)/version.o \
             -L$(CASALIB) -static \
	     -latnf -limages -lms -lcomponents -lcoordinates \
	     -llattices -lfits -lmeasures -lmeasures_f \
	     -ltables -lscimath -lscimath_f -lcasa  \
	     $(CASALIB)/libwcs.a \
	      -lrpfits $(CFITSIOLIB) -dy -llapack -lblas -lg2c

LD        := $(CXX)
LDFLAGS   := -shared -Wl$(TARGET) -s

PYVERSION := 2.3
PYTHONROOT := /usr
PYTHONINC := -I$(PYTHONROOT)/include/python$(PYVERSION)
PYTHONLIB := -L$(PYTHONROOT)/lib -lpython$(PYVERSION) 


BOOSTROOT := /usr
BOOSTLIB  := -L$(BOOSTROOT)/lib -lboost_python
BOOSTINC  := -I$(BOOSTROOT)/include/boost


INCDIRS   := -I. $(CASAINC) $(BOOSTINC) $(PYTHONINC) $(CFITSIOINC)
LIBS      := $(PYTHONLIB)  $(BOOSTLIB) $(CASAPPLIB)

CCOBJECTS := MathUtils.o \
             MathUtils2.o \
	     SDLog.o \
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
	     SDLineFinder.o 
PYOBJECTS := python_SDMemTable.o \
             python_SDReader.o \
             python_SDWriter.o \
             python_SDMath.o \
             python_SDFitter.o \
             python_SDFitTable.o \
	     python_SDLineFinder.o \
	     python_SDLog.o \
             python_SD.o

OBJECTS := $(CCOBJECTS) $(PYOBJECTS)

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
	     SDLog.h \
             python_SD.h

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
