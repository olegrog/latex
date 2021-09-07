PROJECT = $(shell git ls-tree head --name-only | grep tex | cut -f 1 -d '.')
SOURCES = $(PROJECT).tex

%.pdf: $(SOURCES)
	latexmk -xelatex $*

all: $(PROJECT).pdf

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf,snm,nav}

.PHONY: clean
