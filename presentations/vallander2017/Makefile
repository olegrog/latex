DIRS = couette_ couette2_ elliptic_ noncoaxial_
TEX := $(shell ls *.tex)
PDF := $(TEX:.tex=.pdf)

all: $(DIRS) $(PDF)

%_:
	ln -sf ../seminar_cc2016/$* $*

%.pdf: %.tex
	latexmk -pdf -pdflatex='pdflatex --shell-escape' $<

