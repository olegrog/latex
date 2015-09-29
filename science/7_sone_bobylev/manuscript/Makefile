FILE=manuscript
SHELL=/bin/bash
FIG=Fig*.pdf
SOURCES=$(FILE).tex $(FILE).bib $(FILE).rty $(FIG)
CONVERT=gs -sDEVICE=pdfwrite -sProcessColorModel=DeviceGray -sColorConversionStrategy=Gray -dOverrideICC

all: $(FILE).zip $(FILE).pdf

$(FILE).zip:  $(SOURCES)
	zip -r $@ $(SOURCES)

$(FILE).pdf: $(FILE).tex $(FILE).bib $(FIG)
	latex $(FILE)
	bibtex $(FILE)
	latex $(FILE)
	latex $(FILE)

pics:
	cp ../tikz/geometry.pdf             Fig01.pdf~
	cp ../contours/heat-temp.pdf        Fig02.pdf
	cp ../contours/snit-0-temp.pdf      Fig03.pdf
	cp ../contours/nonslip-0-flow.pdf   Fig04.pdf
	cp ../contours/snit-0-flow.pdf      Fig05.pdf
	cp ../contours/asym-0.01-temp.pdf   Fig06.pdf
	cp ../contours/data-0.01-temp.pdf   Fig07.pdf
	cp ../contours/asym-0.01-flow.pdf   Fig08.pdf
	cp ../contours/data-0.01-flow.pdf   Fig09.pdf
	cp ../contours/data-0.1-temp.pdf    Fig10.pdf
	cp ../contours/data-0.1-flow.pdf    Fig11.pdf
	cp ../vs_kn/bottom_temp.pdf         Fig12.pdf
	cp ../vs_kn/bottom_flow.pdf         Fig13.pdf
	cp ../vs_kn/top_temp.pdf            Fig14.pdf
	cp ../vs_kn/top_flow.pdf            Fig15.pdf

Fig%.pdf: Fig%.pdf~
	$(CONVERT) -o $@ -f $<
	rm $<

clean:
	rm -f $(FILE).{aux,bbl,blg,bcf,fff,lof,log,spl,toc,out,pdf,zip}
