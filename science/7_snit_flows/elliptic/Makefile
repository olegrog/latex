FIGURES := $(shell grep 'pdf:' Makefile | sed 's/:.*//' | grep -v FIGURES)

all: $(FIGURES)

asym-0.02-flow.pdf: plot_field.py vtk/asym-0.02.vtk
	./$^ $@ U --factor='(np.sqrt(np.pi/2)*0.02)**-1' --lcmin=1.2,0.3 --lcmax=0.3,0.85 --lpow=0.8

kgf-0.02-flow.pdf: plot_field.py vtk/asym-0.02.vtk
	./$^ $@ U1 --factor='3.42792**-1' --latex='v_{i1}' --lcmin=1.2,0.3 --lcmax=0.3,0.85 --lpow=0.8

kes-0.02-flow.pdf: plot_field.py vtk/kes-0.02.vtk
	./$^ $@ UMean --factor='(np.sqrt(np.pi/2)*0.02)**-1' --kn=0.02 --lcmin=1.3,0.4 --lcmax=0.8,0.8 --lpow=1.2 --lsteps-max=8

clean:
	rm -f $(FIGURES)
