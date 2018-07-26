PROJECT=glasgow
PICS=pics

SOURCES=$(PROJECT).tex

$(PROJECT).pdf: $(SOURCES)
	latex $(PROJECT)
	latex $(PROJECT)

pics:
	mkdir -p elliptic
	cd elliptic && ln -sf ../../../science/7_snit_flows/elliptic/*.pdf .

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf,snm,nav}

.PHONY: pics clean
