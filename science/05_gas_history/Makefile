FILE=history
BIB=first.bib second.bib
TEX=$(FILE).tex titlepage.tex formatting.tex
FIG=mipt_logo.pdf
SOURCES=$(TEX) $(BIB) $(FIG)
SHELL=/bin/bash

all: $(FILE).pdf

$(FILE).pdf: $(SOURCES) 
	latex $(FILE)
	biber $(FILE)
	latex $(FILE)

clean:
	rm -f $(FILE).{aux,lof,toc}
	rm -f $(FILE).{bbl,blg,bcf}
	rm -f $(FILE).fff
	rm -f $(FILE).lof
	rm -f $(FILE).log
	rm -f $(FILE).spl
	rm -f $(FILE).bcf
	rm -f $(FILE).out
	rm -f $(FILE).pdf

