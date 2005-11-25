PREFIX := /usr
PYDIR := $(PREFIX)/lib/python2.3/site-packages

ifndef ASAPROOT
   ASAPROOT := $(shell pwd)
endif

PY := $(wildcard python/*.py)
LIBS := /tmp/_asap.so
BINS := bin/asap

CASAROOT  := $(word 1, $(AIPSPATH))
PLATFORM  := $(word 2, $(AIPSPATH))
DISTDIR   := asap_$(PLATFORM)

all: module #doc

module:
	@cd $(ASAPROOT)/src; make

doc:
	@cd $(ASAPROOT)/doc; make

install:
	@if ( test ! -d $(PYDIR)/asap ) ; then mkdir -p $(PYDIR)/asap ; fi
	@if ( test ! -d $(PREFIX)/bin ) ; then mkdir -p $(PREFIX)/bin ; fi
	@for file in $(LIBS) ; do cp -f $$file $(PYDIR)/asap ; done
	@for file in $(BINS) ; do cp -f $$file $(PREFIX)/bin ; done
	@for file in $(PY) ; do cp -f $$file $(PYDIR)/asap ; done
	@if ( test ! -d $(PREFIX)/share/asap ) ; then mkdir -p $(PREFIX)/share/asap ; fi
	@cp -f share/ipythonrc-asap $(PREFIX)/share/asap/
	@echo "Successfully installed asap module to" $(PYDIR)

clean:
	@cd $(ASAPROOT)/src; make clean
	@cd $(ASAPROOT)/doc; make clean

datadist:
	@echo "Generating ASAP data archive from aips++ installation..."
	@cd $(CASAROOT); tar cfj $(ASAPROOT)/$(DISTDIR)/share/data.tar.bz2 data/ephemerides data/geodetic
	@echo "...done."


dist: module doc
	@cd $(ASAPROOT)
	@if ( test -d $(DISTDIR)  ) ; then rm -rf $(DISTDIR) ; fi
	@mkdir $(DISTDIR)
	@mkdir $(DISTDIR)/build $(DISTDIR)/bin $(DISTDIR)/share
	@for file in $(LIBS) ; do cp -f $$file $(DISTDIR)/build/ ; done
	@for file in $(PY) ; do cp -f $$file $(DISTDIR)/build/ ; done
	@for file in $(BINS) ; do cp -f $$file $(DISTDIR)/bin/ ; done
	@cp -f share/ipythonrc-asap $(DISTDIR)/share/
	make datadist
	@cp -f doc/README $(DISTDIR)/
	@cp -f admin/install.sh $(DISTDIR)/bin/
	@echo "Creating compressed archive..."
	@tar jcf $(DISTDIR).tar.bz2 $(DISTDIR)
	@rm -rf $(DISTDIR)/
	@echo "Successfully created binary package" $(DISTDIR).tar.bz2

.PHONY: install clean
