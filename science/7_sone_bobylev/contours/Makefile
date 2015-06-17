FIGURES := $(shell grep 'pdf:' Makefile | sed 's/:.*//' | grep -v FIGURES)
DATA_DIR = ~/kesolver/data/sone_bobylev
ASYM_DIR = ~/my_openfoam/run/sone_bobylev/snit_solution

all: $(FIGURES)

asym-%.vtu: knudsen_layer_correction.py snit-%.vtk
	./$^ $@ $(*F)

snit-0-temp.pdf: plot_field.py vtk/snit-0.vtk
	./$^ $@ T0 T

snit-0-vel.pdf: plot_field.py vtk/snit-0.vtk
	./$^ $@ U1 v_{i1} 1.148847 -1 1.7

snit-0.01-temp.pdf: plot_field.py vtk/snit-0.01.vtk
	./$^ $@ T0 T_0

asym-0.005-temp.pdf: plot_field.py vtk/asym-0.005.vtu
	./$^ $@ T T

asym-0.01-temp.pdf: plot_field.py vtk/asym-0.01.vtu
	./$^ $@ T T

asym-0.02-temp.pdf: plot_field.py vtk/asym-0.02.vtu
	./$^ $@ T T

asym-0.05-temp.pdf: plot_field.py vtk/asym-0.05.vtu
	./$^ $@ T T

asym-0.01-rho.pdf: plot_field.py vtk/asym-0.01.vtu
	./$^ $@ rho \\rho

asym-0.01-vel.pdf: plot_field.py vtk/asym-0.01.vtu
	./$^ $@ U \\frac{v_i}{k} 0.01 -1

asym-0.05-vel.pdf: plot_field.py vtk/asym-0.05.vtu
	./$^ $@ U \\frac{v_i}{k} 0.05 -1

data-0.01-temp.pdf: plot_field.py vtk/data-0.01.vtk
	./$^ $@ T T 0.01

clean:
	rm -f $(FIGURES)
