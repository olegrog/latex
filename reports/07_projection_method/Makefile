SHELL = /bin/bash
PROJECT = $(patsubst %.tex,%,$(wildcard *.tex))
SOURCES = $(PROJECT).tex

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -pdf $(PROJECT)

.PHONY: clean

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{spl,bbl}
