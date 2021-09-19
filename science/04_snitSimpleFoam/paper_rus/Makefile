PROJECT = paper
PICS = _pics
SOURCES = $(PROJECT).tex $(PICS) info.txt copyright_transfer.doc

all: $(PROJECT).zip $(PROJECT).pdf

$(PROJECT).zip:	$(SOURCES)
	zip -r $@ $(SOURCES)

$(PROJECT).pdf: $(PROJECT).tex $(PICS)
	latexmk -pdf $(PROJECT)

.PHONY: pics clean

pics:
	mkdir -p $(PICS)
	cd $(PICS) && ln -sf ../../transport/Y1_2.pdf Y1.pdf
	cd $(PICS) && ln -sf ../../tikz/geometry.pdf bobylev.pdf
	cd $(PICS) && ln -sf ../../plane/T-asym.pdf T_asym.pdf
	cd $(PICS) && ln -sf ../../plane/T-heat.pdf T_heat.pdf
	cd $(PICS) && ln -sf ../../plane/U.pdf U_0.pdf
	cd $(PICS) && ln -sf ../../plane/U01.pdf U_001.pdf
	cd $(PICS) && ln -sf ../../noncoaxial/cylinder5.pdf cylinders.pdf
	cd $(PICS) && ln -sf ../../noncoaxial/sphere5.pdf spheres.pdf
	cd $(PICS) && ln -sf ../../temper/U-tau-cylinders.pdf alpha.pdf
	cd $(PICS) && ln -sf ../../elliptic/U.pdf elliptic.pdf
	cd $(PICS) && ln -sf ../../elliptic/delta-x.pdf temper.pdf



clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{zip,run.xml,bbl}
