DIR = pics
#FIGURES := $(shell grep 'pdf:' Makefile | sed 's/:.*//' | grep -v FIGURES)
CASES := $(shell ls -d */ | grep '\.' | cut -f1 -d'/')
FIGURES := $(addprefix $(DIR)/, $(addsuffix .pdf, $(CASES)))

all: $(FIGURES)

$(DIR)/%.pdf: plot_field.py %/VTK/*.vtk
	./$^ $@ U --kn=$(lastword $(subst _, ,$*)) --factor='np.sqrt(0.6)' --lrad=.2 --lcmin=-.75,0 --lcmax=-.3,2 --lpow=1.5

$(DIR)/old_1.89080.pdf: plot_field.py old_1.89080/VTK/*.vtk
	./$^ $@ U --kn=1.89080 --factor='np.sqrt(0.6)' --lrad=.5 --lcmin=.5,0 --lcmax=1.2,1.2 --lpow=2

clean:
	rm -f $(FIGURES)
