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
PYRAPDIR := $(EXTDIR)/pyrap-0.3.2/pyrap/Converters

# the casa environment AIPSPATH has to be defined
CASAROOT := $(word 1, $(CASAPATH))
CASAARCH := $(word 2, $(CASAPATH))
ifeq "$(CASAARCH)" "darwin"
   CASADIST := /usr/lib/casapy
else
   CASADIST := /opt/casa
endif


# CASACORE 
COREINCD := $(CASAROOT)/$(CASAARCH)/include/casacore
CORELIBD := $(CASAROOT)/$(CASAARCH)/lib

# boost include directory
BOOSTINCD := 
BOOSTLIBD :=

# numpy include directory (assume to use numpy in CASA)
NUMPYROOT := $(CASADIST)/lib/python2.5/site-packages/numpy
NUMPYINCD := $(NUMPYROOT)/core/include
NUMPYLIBD := $(NUMPYROOT)/core/



OBJECTS := $(PYRAPDIR)/PycArray.o \
           $(PYRAPDIR)/PycArrayNP.o \
           $(PYRAPDIR)/PycBasicData.o \
           $(PYRAPDIR)/PycExp.o \
           $(PYRAPDIR)/PycRecord.o \
           $(PYRAPDIR)/PycValueHolder.o #\
#           $(PYRAPDIR)/PycArrayNA.o

all:
