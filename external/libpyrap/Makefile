#####
# Makefile for libpyrap
#
# ASSUMPTION:
#    CASA+CASACORE build environment
#    numpy is installed (installing CASA is fine since it contains numpy)
#    boost is installed 
#
#####
EXTDIR := $(shell pwd)
PYRAPROOT := $(EXTDIR)/pyrap-0.3.2/ 
PYRAPDIR := $(EXTDIR)/pyrap-0.3.2/pyrap/Converters
TESTDIR := $(PYRAPDIR)/test

# the casa environment AIPSPATH has to be defined
CASAROOT := $(word 1, $(CASAPATH))
CASAARCH := $(word 2, $(CASAPATH))
ifeq "$(CASAARCH)" "darwin"
   CASADIST := /opt/casa
else
   CASADIST := /usr/lib/casapy
endif

TARGET := libpyrap.so
TESTMOD := $(TESTDIR)/_tConvert.so
INSTDIR := $(CASAROOT)/$(CASAARCH)/lib

# CASACORE 
#COREINCD := $(shell $(shell make -s -f $(CASAROOT)/$(CASAARCH)/makedefs VARS=COREINCD eval_vars) && echo $$COREINCD)
COREINCD := $(CASAROOT)/$(CASAARCH)/include/casacore
CORELIBD := $(CASAROOT)/$(CASAARCH)/lib
CORELIB := -L$(CORELIBD) -lcasa_casa

# boost 
BOOSTROOT := $(shell $(shell make -s -f $(CASAROOT)/$(CASAARCH)/makedefs VARS=BOOSTROOT eval_vars) && echo $$BOOSTROOT)
BOOSTINCD := $(BOOSTROOT)/include
BOOSTLIBD := $(BOOSTROOT)/lib
BOOSTLIB := -L$(BOOSTLIBD) -lboost_python

# numpy (use numpy in CASA)
NUMPYROOT := $(CASADIST)/lib/python2.5/site-packages/numpy
NUMPYINCD := $(NUMPYROOT)/core/include
NUMPYLIBD := $(NUMPYROOT)/core/
NUMPYLIB := #-L$(NUMPYLIBD) -lmultiarray 

# python
PYTHONROOT := $(shell $(shell make -s -f $(CASAROOT)/$(CASAARCH)/makedefs VARS=PYTHONROOT eval_vars) && echo $$PYTHONROOT)
PYTHONINCD := $(PYTHONROOT)/include/python2.5

# the compiler
ifndef CXX
   CXX      := g++
endif

# the linker
LD := $(CXX)

# compiler flags
CXXFLAGS := -fPIC -O3 -g
CXXFLAGS += -ansi -Wno-long-long -Wall
CXXOPTS := -DAIPS_64B -DAIPS_USENUMPY

# links to external libraries
LIBS := $(CORELIB) $(BOOSTLIB) $(NUMPYLIB)

# include directories
INCDIRS := -I$(COREINCD) -I$(BOOSTINCD) -I$(NUMPYINCD) \
           -I$(PYRAPROOT) -I$(PYTHONINCD)

# link option
LDFLAGS := -shared 

OBJECTS := $(PYRAPDIR)/PycArray.o \
           $(PYRAPDIR)/PycArrayNP.o \
           $(PYRAPDIR)/PycBasicData.o \
           $(PYRAPDIR)/PycExcp.o \
           $(PYRAPDIR)/PycRecord.o \
           $(PYRAPDIR)/PycValueHolder.o #\
#           $(PYRAPDIR)/PycArrayNA.o

HEADERS := $(PYRAPDIR)/PycArrayComCC.h \
	   $(PYRAPDIR)/PycArrayComH.h \
	   $(PYRAPDIR)/PycArray.h \
	   $(PYRAPDIR)/PycArrayNP.h \
	   $(PYRAPDIR)/PycBasicData.h \
	   $(PYRAPDIR)/PycExcp.h \
	   $(PYRAPDIR)/PycRecord.h \
	   $(PYRAPDIR)/PycValueHolder.h #\
#	   $(PYRAPDIR)/PycArrayNA.h

TESTSRC := $(TESTDIR)/tConvert.cc
TESTOBJ := $(TESTDIR)/tConvert.o

all: $(TARGET)

.cc.o:
	$(CXX) $(CXXOPTS) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

$(TARGET): $(OBJECTS)
	$(LD) $(LDFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

$(OBJECTS): $(HEADERS)

clean:
	rm -rf $(OBJECTS) $(TARGET)
	rm -rf $(TESTDIR)/$(TESTMOD) $(TESTDIR)/$(TESTOBJ)

install:
	cp -f $(TARGET) $(INSTDIR)

test:
	$(CXX) $(CXXOPTS) -c $(CXXFLAGS) $(INCDIRS) -o $(TESTOBJ) $(TESTSRC)
	$(LD) $(LDFLAGS) -o $(TESTMOD) $(TESTOBJ) $(LIBS) -L$(INSTDIR) -lpyrap
