PROJECT = $(shell git ls-files | grep \.tex$$ | cut -f 1 -d '.')
SOURCES = $(PROJECT).tex

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -xelatex $(PROJECT) -shell-escape

.PHONY: clean

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{run.xml,bbl}

