PICS = pics
CASES := $(shell ls -d */ | grep 'km' | cut -f1 -d'/')
FIGURES := $(addprefix $(PICS)/, $(addsuffix .pdf, $(CASES)))

all: $(FIGURES)

$(PICS)/%.pdf: plot_field.py %/VTK/*.vtk
	./$^ $@ U --lrad=3 --factor='np.sqrt(1.2)' --lcmax=0,2

clean:
	rm -f $(FIGURES)
