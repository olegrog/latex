FILE=manuscript
SHELL=/bin/bash
FIG=Fig*.pdf
SOURCES=$(FILE).tex $(FILE).bib $(FILE).rty $(FIG)
CONVERT=gs -sDEVICE=pdfwrite -sProcessColorModel=DeviceGray -sColorConversionStrategy=Gray -dOverrideICC

all: $(FILE).zip $(FILE).pdf

$(FILE).zip:  $(SOURCES)
	zip -r $@ $(SOURCES)

_$(FILE).%: $(FILE).%
	sed "s/$(FILE)/_$(FILE)/" $< | iconv -f UTF-8 -t CP866 > $@

$(FILE).pdf: _$(FILE).tex _$(FILE).bib _$(FILE).rty $(FIG)
	latex _$(FILE)
	bibtex _$(FILE)
	latex _$(FILE)
	latex _$(FILE)
	mv _$@ $@

pics:
	cp ../tikz/sone_bobylev.pdf         Fig01.pdf~
	cp ../contours/heat-temp.pdf        Fig02.pdf
	cp ../contours/snif-0-temp.pdf      Fig03.pdf
	cp ../contours/nonslip-0-flow.pdf   Fig04.pdf
	cp ../contours/snif-0-flow.pdf      Fig05.pdf
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
	cp ../tikz/friedlander.pdf          Fig16.pdf

Fig%.pdf: Fig%.pdf~
	$(CONVERT) -o $@ -f $<
	rm $<

clean:
	rm -f _$(FILE).{aux,bbl,blg,bcf,fff,lof,log,spl,toc,out,tex,bib,pdf}
	rm -f $(FILE).{pdf,zip}
