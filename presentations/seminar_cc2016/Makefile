PROJECT=seminar_mei
PICS=pics

SOURCES=$(PROJECT).tex

$(PROJECT).pdf: $(SOURCES) 
	latex $(PROJECT)
	latex $(PROJECT)

pics:
	mkdir -p couette couette2 snit_flows noncoaxial elliptic
	ln -s ../../science/6_couette_flow/tikz/geometry.pdf couette/
	ln -s ../../science/6_couette_flow/distrib_f/kn*pdf couette/
	ln -s ../../science/6_couette_flow/profiles/*.pdf couette/
	ln -s ../../science/6_couette_flow/integrated/*.pdf couette2/
	ln -s ../../science/7_snit_flows/tikz/sone_bobylev.pdf snit_flows/geometry.pdf
	ln -s ../../science/7_snit_flows/contours/*.pdf snit_flows/
	ln -s ../../science/7_snit_flows/vs_kn/*.pdf snit_flows/
	ln -s ../../science/4_snitSimpleFoam/noncoaxial/*.pdf noncoaxial/
	ln -s ../../science/4_snitSimpleFoam/temper/*.pdf noncoaxial/
	ln -s ../../science/4_snitSimpleFoam/force/*.pdf noncoaxial/
	ln -s ../../science/4_snitSimpleFoam/elliptic/U.pdf elliptic/
	ln -s ../../science/4_snitSimpleFoam/moment/*.pdf elliptic/
	ln -s ../../science/7_snit_flows/elliptic/*.pdf elliptic/
	ln -s ../../science/7_snit_flows/profiles/*.pdf elliptic/

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf,snm,nav}

.PHONY: pics clean
