TARGET   := ../lib/atnf_sd.so

CXX      := g++-3.2.3
CXXFLAGS := -O -fPIC
#CXXFLAGS += -fpermissive # This should be removed.
#CXXFLAGS += -ansi
CXXFLAGS += -pedantic -Wno-long-long
#CXXFLAGS += -Wall -Wno-parentheses

CFITSIOLIB := /usr/local/lib/libcfitsio.a #-lcfitsio
CFITSIOINC := -I/usr/local/include/cfitsio

CASAROOT  := $(word 1, $(AIPSPATH))
CASAINC   := -I$(CASAROOT)/code/include
CASALIB   := $(CASAROOT)/linux_gnu/lib
CASAPPLIB := $(CASALIB)/version.o \
             -L$(CASALIB) \
	     -latnf -lms -limages -lcomponents -lcoordinates \
	     -llattices -lfits -lmeasures -lmeasures_f \
	     -ltables -lscimath -lscimath_f -lcasa \
              $(CASALIB)/libwcs.a -lmir \
             -lrpfits $(CFITSIOLIB) -llapack -lblas -lg2c

LD        := g++-3.2.3
LDFLAGS   := -shared -Wl,--soname=$(TARGET)
LDFLAGS   += -s

BOOSTLIB  := /nfs/atapplic/boost/linux/lib/libboost_python2.3-gcc-mt.a
BOOSTINC  :=

PYTHONINC := -I/usr/local/include/python2.3
PYTHONLIB := -lpython2.3

INCDIRS   := -I. $(CASAINC) $(BOOSTINC) $(PYTHONINC) $(CFITSIOINC)
LIBS      := $(BOOSTLIB) $(PYTHONLIB) $(CASAPPLIB)

OBJECTS   := MathUtils.o \
             SDContainer.o \
             SDMemTable.o \
             SDReader.o \
             SDWriter.o \
             SDMath.o \
             SDTemplates.o \
             python_SDMemTable.o \
             python_SDReader.o \
             python_SDWriter.o \
             python_SDMath.o \
             python_SD.o


HEADERS   := SDContainer.h \
             SDMemTable.h \
             SDReader.h \
             SDWriter.h \
             SDMath.h \
             SDMemTableWrapper.h \
             SDMathWrapper.h \
             SDReaderWrapper.h \
             SDWriterWrapper.h \
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
