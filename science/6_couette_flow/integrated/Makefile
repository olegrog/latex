FIGURES = shear.pdf flow.pdf qflow.pdf qflowy.pdf pxx.pdf pzz.pdf temp.pdf press.pdf diff.pdf
DATA = $(shell ls data-*)
DIFF1 = $(subst data, _diff1, $(DATA))
DIFF2 = $(subst data, _diff2, $(DATA))
DIFF3 = $(subst data, _diff3, $(DATA))

all: $(FIGURES) $(DIFF1) $(DIFF2) $(DIFF3)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: %.sh
	./$<

_diff1-%.txt: asym-%.txt
	./compare.py $* data asym > $@

_diff2-%.txt: dsmc-%.txt
	./compare.py $* dsmc data > $@

_diff3-%.txt: dsmc-%.txt
	./compare.py $* dsmc asym > $@

clean:
	rm -f $(FIGURES)
