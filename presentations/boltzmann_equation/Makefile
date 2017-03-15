PROJECT=boltzmann
PICS=pics

SOURCES=$(PROJECT).tex

$(PROJECT).pdf: $(SOURCES) 
	latexmk -pdf $(PROJECT)

pics:
	mkdir -p elliptic
	cp ../../science/7_snit_flows/elliptic/*.pdf elliptic/

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf,snm,nav,fls,fdb_latexmk}

.PHONY: pics clean
