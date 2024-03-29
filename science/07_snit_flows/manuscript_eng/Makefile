SHELL = /bin/bash
PROJECT = $(shell git ls-tree head --name-only | grep tex$$ | cut -f 1 -d '.')
FIG = Fig*.pdf
SOURCES = $(PROJECT).tex $(PROJECT).bib $(PROJECT).bbl $(FIG)
CONVERT = gs -sDEVICE=pdfwrite -sProcessColorModel=DeviceGray -sColorConversionStrategy=Gray -dOverrideICC

all: $(PROJECT).zip $(PROJECT).pdf

$(PROJECT).zip:  $(SOURCES)
	zip -r $@ $(SOURCES)

$(PROJECT).bbl: $(PROJECT).pdf

%.pdf: %.tex %.bib $(FIG)
	latexmk -pdf $*

.PHONY: pics clean

pics:
	cp ../tikz/sone_bobylev.pdf         Fig01.pdf
	cp ../contours/heat-0-temp.pdf      Fig02.pdf
	cp ../contours/asym-0-temp.pdf      Fig03.pdf
	cp ../contours/nonslip-0-vel.pdf    Fig04.pdf
	cp ../contours/asym-0-vel.pdf       Fig05.pdf
	cp ../contours/asym-0.01-temp.pdf   Fig06.pdf
	cp ../contours/kes-0.01-temp.pdf    Fig07.pdf
	cp ../contours/asym-0.01-vel.pdf    Fig08.pdf
	cp ../contours/kes-0.01-vel.pdf     Fig09.pdf
	cp ../contours/kes-0.1-temp.pdf     Fig10.pdf
	cp ../contours/kes-0.1-vel.pdf      Fig11.pdf
	cp ../vs_kn/bottom_temp.pdf         Fig12.pdf
	cp ../vs_kn/bottom_flow.pdf         Fig13.pdf
	cp ../vs_kn/top_temp.pdf            Fig14.pdf
	cp ../vs_kn/top_flow.pdf            Fig15.pdf
	cp ../elliptic/kgf-0.02-flow.pdf    Fig16.pdf
	cp ../elliptic/second-0.02-flow.pdf Fig17.pdf
	cp ../elliptic/kes-0.02-flow.pdf    Fig18.pdf
	cp ../profiles/outer_temp.pdf       Fig19.pdf
	cp ../profiles/inner_temp.pdf       Fig20.pdf
	cp ../profiles/bottom_temp.pdf      Fig21.pdf
	cp ../profiles/left_temp.pdf        Fig22.pdf
	cp ../profiles/outer_vel.pdf        Fig23.pdf
	cp ../profiles/inner_vel.pdf        Fig24.pdf
	cp ../profiles/bottom_vel.pdf       Fig25.pdf
	cp ../profiles/left_vel.pdf         Fig26.pdf

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{zip,run.xml,bbl}
