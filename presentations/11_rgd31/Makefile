PROJECT=glasgow
PICS=pics

SOURCES=$(PROJECT).tex

$(PROJECT).pdf: $(SOURCES)
	latexmk -pdf -xelatex $(PROJECT)

.PHONY: pics clean

pics:
	mkdir -p elliptic
	cd elliptic && ln -sf ../../../science/07_snit_flows/elliptic/*.pdf .
	ln -sf ../../science/08_vs_gkua/rgd31/cylinder-70km .
	ln -sf ../../science/08_vs_gkua/rgd31/cylinder-90km .

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{bbl}

