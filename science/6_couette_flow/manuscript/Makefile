FILE=manuscript
SHELL=/bin/bash
FIG=Fig*.pdf
SOURCES=$(FILE).tex $(FILE).bib highlights.txt $(FIG)

all: $(FILE).zip $(FILE).pdf

$(FILE).zip:  $(SOURCES)
	zip -r $@ $(SOURCES)

$(FILE).pdf: $(FILE).tex $(FILE).bib $(FIG)
	latex $(FILE)
	bibtex $(FILE)
	latex $(FILE)
	latex $(FILE)

pics:
	#cp ../exclusions/max.pdf                Fig01.pdf
	#cp ../exclusions/mean.pdf               Fig02.pdf
	cp ../distrib_f/kn0.1-boundary.pdf      Fig03.pdf
	cp ../distrib_f/kn0.1-center.pdf        Fig04.pdf
	cp ../distrib_f/kn1.0-boundary.pdf      Fig05.pdf
	cp ../distrib_f/kn1.0-center.pdf        Fig06.pdf
	cp ../distrib_f/kn10-boundary.pdf       Fig07.pdf
	cp ../distrib_f/kn10-center.pdf         Fig08.pdf
	cp ../profiles/Pxy.pdf                  Fig09.pdf
	cp ../profiles/vx.pdf                   Fig10.pdf
	cp ../profiles/qx.pdf                   Fig11.pdf
	cp ../profiles/qy.pdf                   Fig12.pdf
	cp ../profiles/P.pdf                    Fig13.pdf
	cp ../profiles/Pyy.pdf                  Fig14.pdf
	cp ../profiles/Pxx.pdf                  Fig15.pdf
	cp ../profiles/Pzz.pdf                  Fig16.pdf
	cp ../profiles/tau.pdf                  Fig17.pdf
	cp ../profiles/omega.pdf                Fig18.pdf
	cp ../integrated/shear.pdf              Fig19.pdf
	cp ../integrated/flow.pdf               Fig20.pdf
	cp ../integrated/qflow.pdf              Fig21.pdf
	cp ../integrated/qflowy.pdf             Fig22.pdf
	cp ../integrated/pxx.pdf                Fig23.pdf
	cp ../integrated/pzz.pdf                Fig24.pdf
	cp ../integrated/temp.pdf               Fig25.pdf
	cp ../integrated/diff.pdf               Fig26.pdf
	#cp ../gamma/B.pdf                       Fig27.pdf
	cp ../gamma/B_4.pdf                     Fig28.pdf
	cp ../gamma/T1_1.pdf                    Fig29.pdf
	cp ../gamma/T1_2.pdf                    Fig30.pdf
	cp ../gamma/T2_1.pdf                    Fig31.pdf
	cp ../gamma/T2_2.pdf                    Fig32.pdf
	cp ../gamma/TT_12.pdf                   Fig33.pdf
	cp ../gamma/TT_2.pdf                    Fig34.pdf
	cp ../gamma/Q_2.pdf                     Fig35.pdf
	cp ../gamma/Q_3.pdf                     Fig36.pdf
	cp ../gamma/QQ_22.pdf                   Fig37.pdf
	cp ../gamma/QQ_3.pdf                    Fig38.pdf

clean:
	rm -f $(FILE).{aux,bbl,blg,fff,lof,log,spl,toc,pdf,zip}
