PROJECT=boltzmann
PICS=pics

SOURCES=$(PROJECT).tex

$(PROJECT).pdf: $(SOURCES) 
	latexmk -pdfxe $(PROJECT)

pics:
	mkdir -p elliptic
	cp ../../science/07_snit_flows/elliptic/*.pdf elliptic/
	latexmk -pdf layers.tex
	latexmk -pdf spectra.tex

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf,snm,nav,fls,fdb_latexmk}

.PHONY: pics clean
