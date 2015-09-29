FIGURES := $(shell grep 'pdf:' Makefile | sed 's/:.*//' | grep -v FIGURES)
DATA_DIR = ~/kesolver/data/sone_bobylev
ASYM_DIR = ~/my_openfoam/run/sone_bobylev/snit_solution

all: $(FIGURES)

asym-%.vtk: knudsen_layer_correction.py snit-%.vtk
	./$^ $@ $(*F)


heat-temp.pdf: plot_field.py vtk/heat.vtk
	./$^ $@ T0 T_0

nonslip-0-temp.pdf: plot_field.py vtk/nonslip-0.vtk
	./$^ $@ T0 T_0

slip-0-temp.pdf: plot_field.py vtk/slip-0.vtk
	./$^ $@ T0 T_0

snit-0-temp.pdf: plot_field.py vtk/snit-0.vtk
	./$^ $@ T0 T_0

snit-0.01-temp.pdf: plot_field.py vtk/snit-0.01.vtk
	./$^ $@ T0 T_0

asym-0.005-temp.pdf: plot_field.py vtk/asym-0.005.vtk
	./$^ $@ T T_0

asym-0.01-temp.pdf: plot_field.py vtk/asym-0.01.vtk
	./$^ $@ T T_0

asym-0.02-temp.pdf: plot_field.py vtk/asym-0.02.vtk
	./$^ $@ T T_0 --lmax=0.15 --lmin=-0.1

asym-0.05-temp.pdf: plot_field.py vtk/asym-0.05.vtk
	./$^ $@ T T_0 --lmax=0

data-0.01-temp.pdf: plot_field.py vtk/data-0.01.vtk
	./$^ $@ TMean T --kn 0.01

data-0.1-temp.pdf: plot_field.py vtk/data-0.1.vtk
	./$^ $@ TMean T --kn 0.1


asym-0.01-rho.pdf: plot_field.py vtk/asym-0.01.vtk
	./$^ $@ rho \\rho


snit-0-flow.pdf: plot_field.py vtk/snit-0.vtk
	./$^ $@ U1 v_{i1} --factor='1.01814**-1' --lpow=1.7

slip-0-flow.pdf: plot_field.py vtk/slip-0.vtk
	./$^ $@ U1 v_{i1} --factor='1.01814**-1' --lpow=1.7

nonslip-0-flow.pdf: plot_field.py vtk/nonslip-0.vtk
	./$^ $@ U1 v_{i1} --factor='1.00494**-1' --lpow=0.8

asym-0.01-flow.pdf: plot_field.py vtk/asym-0.01.vtk
	./$^ $@ U/k v_{i1}

asym-0.05-flow.pdf: plot_field.py vtk/asym-0.05.vtk
	./$^ $@ U/k v_{i1}

data-0.01-flow.pdf: plot_field.py vtk/data-0.01.vtk
	./$^ $@ UMean \\frac{v_i}{k} --factor='(np.sqrt(np.pi/2)*0.01)**-1' --kn=0.01 --lpow=1.2

data-0.1-flow.pdf: plot_field.py vtk/data-0.1.vtk
	./$^ $@ UMean \\frac{v_i}{k} --factor='(np.sqrt(np.pi/2)*0.01)**-1' --kn=0.1 --lpow=0.9

clean:
	rm -f $(FIGURES)
