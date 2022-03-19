PROJECT = $(shell git ls-tree head --name-only | grep \.tex$$ | cut -f 1 -d '.')
SOURCES = $(PROJECT).tex
SUBDIRS = plots

all: $(PROJECT).pdf

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latexmk -xelatex $(PROJECT) -shell-escape

.PHONY: all clean $(SUBDIRS)

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{run.xml,bbl,gnuploterrors}
	rm -rf gnuplottex

