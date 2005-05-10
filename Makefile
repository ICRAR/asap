PREFIX := /usr/local
PYDIR := $(PREFIX)/lib/python2.3/site-packages

ifndef ASAPROOT
   ASAPROOT := $(shell pwd)
endif

PY := $(wildcard python/*.py)
LIBS := lib/_asap.so
BINS := bin/asap

all: module #doc

module:
	@cd $(ASAPROOT)/src; make

doc:
	@cd $(ASAPROOT)/doc; make

install: 
	@if ( test ! -d $(PYDIR)/asap ) ; then mkdir -p $(PYDIR)/asap ; fi
	@for file in $(LIBS) ; do cp -f $$file $(PYDIR)/asap ; done
	@for file in $(BINS) ; do cp -f $$file $(PREFIX)/bin ; done
	@for file in $(PY) ; do cp -f $$file $(PYDIR)/asap ; done
	@echo "Successfully installed asap module to" $(PYDIR)

clean:
	@cd $(ASAPROOT)/src; make clean
	@cd $(ASAPROOT)/doc; make clean    

.PHONY: install clean
