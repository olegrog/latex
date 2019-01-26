FILE=manuscript
FIG= #Fig*.pdf
SOURCES=$(FILE).tex $(FILE).bib $(FIG)

all: $(FILE).zip $(FILE).pdf

$(FILE).zip:  $(SOURCES)
	zip -r $@ $(SOURCES)

$(FILE).pdf: $(FILE).tex $(FILE).bib $(FIG)
	latexmk -pdf $(FILE)

pics:
	cp ../tikz/sone_bobylev.pdf         Fig01.pdf

clean:
	latexmk -C
	rm -f $(FILE).{pdf,zip}
