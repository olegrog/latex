VTK := $(shell ls vtk | grep asym) $(shell ls vtk | grep kes) $(shell ls vtk | grep nonslip)
FIGURES := $(VTK:.vtk=-temp.pdf) $(VTK:.vtk=-vel.pdf) heat-0-temp.pdf

all: $(FIGURES)

%-0-temp.pdf: plot_field.py vtk/%-0.vtk
	./$^ $@ T0 T_0

asym-0-vel.pdf: plot_field.py vtk/asym-0.vtk
	./$^ $@ U1 v_{i1} --factor='1.018270**-1' --lpow=1.7

nonslip-0-vel.pdf: plot_field.py vtk/nonslip-0.vtk
	./$^ $@ U1 v_{i1} --factor='1.00494**-1' --lpow=0.8

asym-%-temp.pdf: plot_field.py vtk/asym-%.vtk
	./$^ $@ T T

asym-0.02-temp.pdf: plot_field.py vtk/asym-0.02.vtk
	./$^ $@ T T --lmax=0.15 --lmin=-0.1

asym-0.05-temp.pdf: plot_field.py vtk/asym-0.05.vtk
	./$^ $@ T T --lmax=0

kes-%-temp.pdf: plot_field.py vtk/kes-%.vtk
	./$^ $@ TMean T --kn $*

asym-%-vel.pdf: plot_field.py vtk/asym-%.vtk
	./$^ $@ U/k v_{i1}

kes-%-vel.pdf: plot_field.py vtk/kes-%.vtk
	./$^ $@ UMean \\frac{v_i}{k} --factor='(np.sqrt(np.pi/2)*$*)**-1' --kn=$*

asym-%-rho.pdf: plot_field.py vtk/asym-%.vtk
	./$^ $@ rho \\rho

clean:
	rm -f $(FIGURES)
