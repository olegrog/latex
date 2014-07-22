FILE=elsarticle
FIG=Fig*.pdf
SOURCES=$(FILE).tex $(FILE).bib $(FIG) highlights.txt

all: $(FILE).zip $(FILE).pdf

$(FILE).zip:  $(SOURCES)
	zip -r $@ $(SOURCES)

$(FILE).pdf: $(FILE).tex $(FILE).bib $(FIG)
	latex $(FILE)
	bibtex $(FILE)
	latex $(FILE)
	latex $(FILE)
