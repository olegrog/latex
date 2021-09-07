PROJECT = manuscript
SOURCES = $(PROJECT).tex

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latex $(PROJECT)
	latex $(PROJECT)

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf}

