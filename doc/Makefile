TARGET := userguide
WEBLOC := /tmp

LATEXHTML := /usr/bin/latex2html
LATEXHTMLOPTS := -local_icons
DVIPDF := /usr/bin/dvipdf
LATEX := /usr/bin/latex

all: pdf html

$(TARGET): pdf html

latex:
	$(LATEX) $(TARGET)
	makeindex $(TARGET)
	$(LATEX) $(TARGET)

pdf: latex
	$(DVIPDF) $(TARGET)

html:
	$(LATEXHTML) $(LATEXHTMLOPTS) $(TARGET)

clean:
	@rm -f $(TARGET).log $(TARGET).aux $(TARGET).dvi

realclean: clean
	@rm -f $(TARGET).pdf $(TARGET)

install:
	-cp -rf  $(TARGET) $(WEBLOC)
	-cp -f $(TARGET).pdf $(WEBLOC)

.PHONY: install clean realclean
