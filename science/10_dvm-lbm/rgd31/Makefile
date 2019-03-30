FILE=manuscript

all: $(FILE).pdf

$(FILE).pdf: $(FILE).tex ../dvm-lbm.bib
	latexmk -pdf $(FILE)

clean:
	latexmk -C
