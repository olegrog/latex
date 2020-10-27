PROJECT = $(shell git ls-tree head --name-only | grep \.tex$$ | cut -f 1 -d '.')
SOURCES = $(PROJECT).tex

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -pdf -xelatex $(PROJECT) -shell-escape

.PHONY: clean

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{run.xml,bbl}
	rm -rf _minted-$(PROJECT)
