TARGET := ./cookbook
WEBLOC := /tmp

LATEXHTML := /usr/bin/latex2html
LATEXHTMLOPTS := -local_icons
DVIPDF := /usr/bin/dvipdf
LATEX := /usr/bin/latex

all: pdf html

$(TARGET): pdf html

latex:
	$(LATEX) $(TARGET).tex

pdf: latex
	$(DVIPDF) $(TARGET).dvi

html:
	$(LATEXHTML) $(LATEXHTMLOPTS) $(TARGET).tex

clean:
	@rm -f $(TARGET).log $(TARGET).aux $(TARGET).dvi

realclean: clean
	@rm -f $(TARGET).pdf $(TARGET)

install:
	-cp -rf  $(TARGET) $(WEBLOC)
	-cp -f $(TARGET).pdf $(WEBLOC)

.PHONY: install clean realclean
