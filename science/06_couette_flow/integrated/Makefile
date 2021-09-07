FIGURES = shear.pdf flow.pdf qflow.pdf qflowy.pdf pxx.pdf pzz.pdf temp.pdf press.pdf diff.pdf
DATA = $(shell ls data-*.txt)
ASYM = $(shell ls asym-*.txt)
DSMC = $(shell ls dsmc-*.txt)
DIFF1 = $(subst data, _diff1, $(DATA))
DIFF2 = $(subst data, _diff2, $(DATA))
DIFF3 = $(subst data, _diff3, $(DATA))

all: $(FIGURES) $(DIFF1) $(DIFF2) $(DIFF3)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: %.sh $(DATA) $(ASYM) $(DSMC)
	./$<

_diff1-%.txt: compare.py asym-%.txt
	./$< $* data asym > $@

_diff2-%.txt: compare.py dsmc-%.txt
	./$< $* dsmc data > $@

_diff3-%.txt: compare.py dsmc-%.txt
	./$< $* dsmc asym > $@

asym-%.txt: attach-knudsen-layer.py ns-%.txt
	./$^ $@

clean:
	rm -f $(FIGURES)
