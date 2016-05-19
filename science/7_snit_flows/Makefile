SHELL = /bin/bash
PROJECT = snit_flows
SUBDIRS = contours vs_kn vs_kn2 tikz
SOURCES = $(PROJECT).tex $(PROJECT).bib

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latex $(PROJECT)
	biber $(PROJECT)
	latex $(PROJECT)

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf}

gsb:
	./discrete_error.py 1.5 12 15 4.5 0.05 0.005

gell:
	./discrete_error.py 5 16 16 8 0.1 0.1
	

