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

pics:
	cp ../transport/Y1_2.pdf            Fig01.pdf
	cp ../tikz/bobylev.pdf              Fig02.pdf
	cp ../moving_wall/T.pdf             Fig03.pdf
	cp ../moving_wall/T-heat.pdf        Fig04.pdf
	cp ../moving_wall/U.pdf             Fig05.pdf
	cp ../moving_wall/U01.pdf           Fig06.pdf
	cp ../noncoaxial/cylinder5.pdf      Fig07.pdf
	cp ../noncoaxial/sphere5.pdf        Fig08.pdf
	cp ../force/forces.pdf              Fig09.pdf
	cp ../force/forces-inv.pdf          Fig10.pdf
	cp ../force/forces-alpha.pdf        Fig11.pdf
	cp ../force/profiles.pdf            Fig12.pdf
	cp ../force/terms-cylinder.pdf      Fig13.pdf
	cp ../force/terms-cylinder-inv.pdf  Fig14.pdf
	cp ../elliptic/U.pdf                Fig15.pdf
	cp ../moment/moment_beta.pdf        Fig16.pdf
	cp ../moment/profiles.pdf           Fig17.pdf
	cp ../elliptic/alpha.pdf            Fig18.pdf
	cp ../elliptic/temper.pdf           Fig19.pdf

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

