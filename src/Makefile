TARGET    = ../lib/atnf_sd.so

CXX	  = g++-3.2.3
CXXFLAGS  = -fPIC -w -g -fpermissive -O2

AIPSROOT = $(word 1, $(AIPSPATH))

AIPSINC = -I$(AIPSROOT)/code/include \
          -I$(AIPSROOT)/code/atnf \
          -I$(AIPSROOT)/code/trial \
          -I$(AIPSROOT)/code/aips

AIPSLIB  = $(AIPSROOT)/linux_gnu/lib
AIPSPPLIB = -L$(AIPSLIB) \
        $(AIPSLIB)/version.o \
        $(AIPSLIB)/libatnf.a \
        $(AIPSLIB)/libtrial.a \
        $(AIPSLIB)/libaips_f.a \
        $(AIPSLIB)/libaips.a \
        -llapack -lblas -lg2c -lcfitsio -lrpfits \
        -lglish -ledit -lsos -lnpd -lwcs

LINK      = g++-3.2.3
LFLAGS    = -shared -Wl,-soname,$(TARGET)
#-rpath,/home/mar637/lib
#LFLAGS    += -Xlinker -rpath -Xlinker -L$(AIPSLIB)

BOOSTLIB  = -lboost_python-mt
BOOSTINC  = 

PYTHONINC = -I/usr/include/python2.2
PYTHONLIB = -lpython2.2

CFITSIOINC = -I/usr/local/include/cfitsio

INCDIRS   = -I. $(AIPSINC) $(BOOSTINC) $(PYTHONINC) $(CFITSIOINC)
LIBDIRS   = $(BOOSTLIB) $(PYTHONLIB) $(AIPSPPLIB)

OBJECTS   = MathUtils.o \
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

HEADERS   = SDContainer.h \
            SDMemTable.h \
            SDReader.h \
            SDWriter.h \
            SDMath.h \
            SDMemTableWrapper.h \
            python_SD.h

.SUFFIXES: .c .o .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cxx.o:	
	$(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

$(TARGET):  $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS)  $(LIBDIRS)

$(OBJECTS): $(HEADERS)

clean:
	-rm -f $(OBJECTS)

realclean: clean
	-rm -f $(TARGET)
