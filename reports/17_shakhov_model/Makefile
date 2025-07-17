PROJECT = $(shell git ls-files | grep \.tex$$ | cut -f 1 -d '.')
SOLVER = relaxation/solver.py
CASES = _piecewise _piecewise_hole _pmaxwell _pmaxwell_Q2 _grad13 _grad13_Q1C _grad13_N4
FIGURES = $(addsuffix /solution.pdf, $(CASES))
SOURCES = $(PROJECT).tex $(PROJECT).bib $(FIGURES)

all: $(PROJECT).pdf $(PROJECT).zip

$(PROJECT).pdf: $(SOURCES)
	latexmk -xelatex $(PROJECT) -shell-escape

$(PROJECT).zip: $(SOURCES)
	zip -r $@ $(PROJECT).tex $(PROJECT).bib $(CASES)

options_piecewise = \
	--function=piecewise --inner=0 --outer=2 --ratio=2 --radius=4
options_piecewise_hole = \
	--function=piecewise --inner=1 --outer=2 --ratio=2 --radius=4
options_pmaxwell = \
	--function=pmaxw --qflow=0.1 --radius=3
options_pmaxwell_Q2 = \
	--function=pmaxw --qflow=2 --radius=5
options_grad13 = \
	--function=grad13 --qflow=0.1 --radius=3
options_grad13_Q1C = \
	--function=grad13,constant --qflow=1 --value=-2e-3 --radius=4
options_grad13_N4 = \
	--function=grad13 --qflow=0.1 --radius=4 -N=4

_%/solution.pdf: $(SOLVER)
	mkdir -p _$*
	cd _$* && ../$(SOLVER) --pdf -N=1000 --end=4.1 --timestep=0.005 --plot=200 --grid=uniform \
		--convolution=full $(options_$*)

.PHONY: clean

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{run.xml,bbl}
	rm -rf *gnuplot*
	rm -rf $(CASES)

