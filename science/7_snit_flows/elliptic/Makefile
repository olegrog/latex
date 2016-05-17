FIGURES := $(shell grep 'pdf:' Makefile | sed 's/:.*//' | grep -v FIGURES)

all: $(FIGURES)

asym-jumps-0.02-flow.pdf: plot_field.py vtk/asym-jumps-0.02.vtk
	./$^ $@ U --factor='(np.sqrt(np.pi)/2*0.02)**-1' --lcmin=1.2,0.2 --lcmax=0.45,0.75 --lpow=0.8

asym-second-0.02-flow.pdf: plot_field.py vtk/asym-second-0.02.vtk
	./$^ $@ U --factor='(np.sqrt(np.pi)/2*0.02)**-1' --lcmin=1.2,0.2 --lcmax=0.45,0.75 --lpow=0.8

kgf-0.02-flow.pdf: plot_field.py vtk/asym-jumps-0.02.vtk
	./$^ $@ U1 --latex='v_{i1}' --lcmin=1.2,0.2 --lcmax=0.35,0.85 --lpow=0.8

kes-0.02-flow.pdf: plot_field.py vtk/kes-0.02.vtk
	./$^ $@ UMean --factor='(np.sqrt(np.pi/2)*0.02)**-1' --kn=0.02 --lcmin=1.3,0.4 --lcmax=0.8,0.8 --lpow=1.2 --lsteps-max=8

clean:
	rm -f $(FIGURES)
