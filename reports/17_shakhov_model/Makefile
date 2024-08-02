PROJECT = $(shell git ls-files | grep \.tex$$ | cut -f 1 -d '.')
SOLVER = relaxation/solver.py
CASES = _piecewise _piecewise_hole _grad13
FIGURES = $(addsuffix /solution.pdf, $(CASES))
SOURCES = $(PROJECT).tex $(FIGURES)

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -xelatex $(PROJECT) -shell-escape

options_piecewise = \
	--function=piecewise --inner=0 --outer=2 --ratio=2
options_piecewise_hole = \
	--function=piecewise --inner=1 --outer=2 --ratio=2
options_grad13 = \
	--function=grad13 --qflow=1 --radius=3

_%/solution.pdf: $(SOLVER)
	mkdir -p _$*
	cd _$* && ../$(SOLVER) --pdf -N=200 --end=5 --timestep=0.05 --plot=20 --grid=uniform $(options_$*)

.PHONY: clean

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{run.xml,bbl}
	rm -rf $(CASES)

