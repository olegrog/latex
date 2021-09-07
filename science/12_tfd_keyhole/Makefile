SHELL = /bin/bash
PROJECT = $(shell git ls-tree head --name-only | grep tex$$ | cut -f 1 -d '.')
SOURCES = $(PROJECT).tex $(PROJECT).bib cases.txt
PICS = _pics

all: pics $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -pdf $(PROJECT) --shell-escape

.PHONY: clean pics zip

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{run.xml,spl,bbl,gnuploterrors}
	rm -rf gnuplottex

pics:
	mkdir -p $(PICS)
	cd $(PICS) && ln -sf ../logos/*.pdf .

zip: $(SOURCES)
	zip -r $(PROJECT).zip $(SOURCES) $(PROJECT).bbl $(PICS)
