
# directoty to install ASAP and the source directory
ifdef NARRABRI_ASAP
ASAPDIR := /DATA/KAPUTAR_2/vor010/ASAP/site-packages/asap
else
ASAPDIR := /usr/local/lib/python2.3/site-packages/asap
endif

SRCDIR := $(shell pwd)

all:
	cd src; make

install: all
	sh - bin/install.sh $(ASAPDIR) $(SRCDIR)
