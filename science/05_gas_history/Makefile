PROJECT = history
BIB = first.bib second.bib
TEX = $(PROJECT).tex titlepage.tex formatting.tex
FIG = mipt_logo.pdf
SOURCES = $(TEX) $(BIB) $(FIG)
SHELL = /bin/bash

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -pdf $(PROJECT)

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{bbl}
