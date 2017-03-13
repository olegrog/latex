PROJECT=seminar_mei
PICS=pics

SOURCES=$(PROJECT).tex

$(PROJECT).pdf: $(SOURCES) 
	latex $(PROJECT)
	latex $(PROJECT)

pics:
	mkdir -p couette couette2 snit_flows noncoaxial elliptic
	cp ../../science/6_couette_flow/tikz/geometry.pdf couette/
	cp ../../science/6_couette_flow/distrib_f/kn*pdf couette/
	cp ../../science/6_couette_flow/profiles/*.pdf couette/
	cp ../../science/6_couette_flow/integrated/*.pdf couette2/
	cp ../../science/7_snit_flows/tikz/snit_flows.pdf snit_flows/geometry.pdf
	cp ../../science/7_snit_flows/contours/*.pdf snit_flows/
	cp ../../science/7_snit_flows/vs_kn/*.pdf snit_flows/
	cp ../../science/4_snitSimpleFoam/noncoaxial/*.pdf noncoaxial/
	cp ../../science/4_snitSimpleFoam/temper/*.pdf noncoaxial/
	cp ../../science/4_snitSimpleFoam/force/*.pdf noncoaxial/
	cp ../../science/4_snitSimpleFoam/elliptic/U.pdf elliptic/
	cp ../../science/4_snitSimpleFoam/moment/*.pdf elliptic/
	cp ../../science/7_snit_flows/elliptic/*.pdf elliptic/

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf,snm,nav}

.PHONY: pics clean
