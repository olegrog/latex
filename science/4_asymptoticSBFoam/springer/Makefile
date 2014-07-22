FILE=springer
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


clean:
	rm -f $(FILE).aux
	rm -f $(FILE).bbl
	rm -f $(FILE).blg
	rm -f $(FILE).fff
	rm -f $(FILE).lof
	rm -f $(FILE).log
	rm -f $(FILE).spl
	rm -f $(FILE).pdf
	rm -f $(FILE).zip

