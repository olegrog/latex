FIGURES := $(shell grep 'pdf:' Makefile | sed 's/:.*//' | grep -v FIGURES)

all: $(FIGURES)

kgf-0.005-flow.pdf: plot_field.py vtk/kgf-0.005.vtk
	./$^ $@ U --factor='(np.sqrt(np.pi/2)*0.005)**-1'

kes-0.005-flow.pdf: plot_field.py vtk/kes-0.005.vtk
	./$^ $@ UMean --factor='(np.sqrt(np.pi/2)*0.005)**-1' --kn=0.005

clean:
	rm -f $(FIGURES)
