FILE=manuscript
FIG=Fig*.pdf
SOURCES=$(FILE).tex $(FILE).bib $(FIG)

all: $(FILE).zip $(FILE).pdf

$(FILE).zip:  $(SOURCES)
	zip -r $@ $(SOURCES)

$(FILE).pdf: $(FILE).tex $(FILE).bib $(FIG)
	latex $(FILE)
	bibtex $(FILE)
	latex $(FILE)
	latex $(FILE)

pics:
	cp ../transport/Y1_2.pdf                    Fig01.pdf
	cp ../tikz/geometry.pdf                     Fig02.pdf
	cp ../plane/T-asym.pdf                      Fig03.pdf
	cp ../plane/T-heat.pdf                      Fig04.pdf
	cp ../plane/U.pdf                           Fig05.pdf
	cp ../plane/U01.pdf                         Fig06.pdf
	cp ../noncoaxial/cylinder5.pdf              Fig07.pdf
	cp ../noncoaxial/sphere5.pdf                Fig08.pdf
	cp ../temper/U-tau-cylinders.pdf            Fig09.pdf
	cp ../force/terms-cylinder-inner.pdf	    Fig10.pdf
	cp ../force/terms-cylinder-outer.pdf        Fig11.pdf
	cp ../force/terms-cylinder-inner-inv.pdf    Fig12.pdf
	cp ../temper/F-tau-cylinders-inner.pdf      Fig13.pdf
	cp ../temper/F-tau-cylinders-outer.pdf      Fig14.pdf
	cp ../force/forces.pdf                      Fig15.pdf
	cp ../force/forces-close.pdf                Fig16.pdf
	cp ../elliptic/U.pdf                        Fig17.pdf
	cp ../moment/moment-beta.pdf                Fig18.pdf
	cp ../moment/profiles.pdf                   Fig19.pdf

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

