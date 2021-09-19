PROJECT=ghost_effect
FIGURES=transport/Y1.pdf fluid/T_asym.pdf fluid/T_heat.pdf fluid/U_asym.pdf fluid/U_kn001.pdf kinetic/T_kn001.pdf kinetic/U_kn001.pdf temper/temper.pdf
SOURCES=$(PROJECT).tex $(PROJECT).bib $(FIGURES)
SHELL=/bin/bash

all: $(PROJECT).pdf

transport/Y1.pdf: transport/Y1.py transport/Y1_2.txt
	cd transport && ./Y1.py

fluid/T_asym.pdf: fluid/temperature.py fluid/asymptotic_N50.txt
	cd fluid && ./temperature.py asymptotic_N50.txt T_asym.pdf

fluid/T_heat.pdf: fluid/temperature.py fluid/heat.txt
	cd fluid && ./temperature.py heat.txt T_heat.pdf

fluid/U_asym.pdf: fluid/velocity.py fluid/asymptotic_N50.txt fluid/asymptotic_N100.txt
	cd fluid && ./velocity.py 0 U_asym.pdf

fluid/U_kn001.pdf: fluid/velocity.py fluid/asymptotic_N50.txt fluid/asymptotic_N100.txt
	cd fluid && ./velocity.py 0.01 U_kn001.pdf

kinetic/T_kn001.pdf: kinetic/temperature.py kinetic/kn001.txt
	cd kinetic && ./temperature.py kn001.txt T_kn001.pdf

kinetic/U_kn001.pdf: kinetic/velocity.py kinetic/kn001.txt
	cd kinetic && ./velocity.py 0.01 kn001.txt U_kn001.pdf

temper/temper.pdf: temper/plot.py temper/1.txt temper/2.txt
	cd temper && ./plot.py temper.pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -pdf $(PROJECT)

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{run.xml,bbl}

clean-pics:
	rm -f {fluid,kinetic,temper}/*.pdf
