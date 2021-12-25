SHELL = /bin/bash
PROJECT = $(patsubst %.tex,%,$(wildcard *.tex))
SOURCES = $(PROJECT).tex

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -xelatex $(PROJECT)

.PHONY: clean

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{bbl}

