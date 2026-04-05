PROJECT=seminar_mei
PICS=pics

SOURCES=$(PROJECT).tex

$(PROJECT).pdf: $(SOURCES) tikz/layers.pdf pics
	latexmk -pdf $(PROJECT)

tikz/layers.pdf: tikz/layers.tex
	cd tikz && latexmk -pdf layers

pics:
	mkdir -p couette couette2 snit_flows noncoaxial elliptic
	ln -sf `realpath ../../science/06_couette_flow/tikz/geometry.pdf` couette/
	ln -sf `realpath ../../science/06_couette_flow/distrib_f/kn*pdf` couette/
	ln -sf `realpath ../../science/06_couette_flow/profiles/*.pdf` couette/
	ln -sf `realpath ../../science/06_couette_flow/integrated/*.pdf` couette2/
	ln -sf `realpath ../../science/07_snit_flows/tikz/sone_bobylev.pdf` snit_flows/geometry.pdf
	ln -sf `realpath ../../science/07_snit_flows/contours/*.pdf` snit_flows/
	ln -sf `realpath ../../science/07_snit_flows/vs_kn/*.pdf` snit_flows/
	ln -sf `realpath ../../science/04_snitSimpleFoam/noncoaxial/*.pdf` noncoaxial/
	ln -sf `realpath ../../science/04_snitSimpleFoam/temper/*.pdf` noncoaxial/
	ln -sf `realpath ../../science/04_snitSimpleFoam/force/*.pdf` noncoaxial/
	ln -sf `realpath ../../science/04_snitSimpleFoam/elliptic/U.pdf` elliptic/
	ln -sf `realpath ../../science/04_snitSimpleFoam/moment/*.pdf` elliptic/
	ln -sf `realpath ../../science/07_snit_flows/elliptic/*.pdf` elliptic/
	ln -sf `realpath ../../science/07_snit_flows/profiles/*.pdf` elliptic/

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{snm,nav}
	rm -r couette/ couette2/ snit_flows/ noncoaxial/ elliptic/

.PHONY: pics clean
