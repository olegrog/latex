FIGURES := $(shell grep 'pdf:' Makefile | sed 's/:.*//' | grep -v FIGURES)

all: $(FIGURES)

asym-0.02-flow.pdf: plot_field.py vtk/kgf-0.02.vtk
	./$^ $@ U --factor='(np.sqrt(np.pi/2)*0.02)**-1' --lpow=0.5

kgf-0.02-flow.pdf: plot_field.py vtk/kgf-0.02.vtk
	./$^ $@ U1 --factor='3.42792**-1' --lpow=0.5

kes-0.02-flow.pdf: plot_field.py vtk/kes-0.02.vtk
	./$^ $@ U --factor='(np.sqrt(np.pi/2)*0.02)**-1' --kn=0.02

clean:
	rm -f $(FIGURES)
