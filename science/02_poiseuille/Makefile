PROJECT=poiseuille
FIGURES=pics/graph.pdf
SOURCES=$(PROJECT).tex $(PROJECT).bib $(FIGURES)
SHELL=/bin/bash
SUBDIRS=pics

all: $(PROJECT).pdf

pics/graph.pdf:	pics/*txt pics/graph.sh
	cd pics && ./graph.sh && pdflatex graph

$(PROJECT).pdf: $(SOURCES) 
	latex $(PROJECT)
	biber $(PROJECT)
	latex $(PROJECT)

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf}
	rm -f pics/graph.{tex,pdf,aux,log} pics/graph-inc*

