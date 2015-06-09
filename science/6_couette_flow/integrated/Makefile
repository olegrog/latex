FIGURES = shear.pdf flow.pdf qflow.pdf qflowy.pdf pxx.pdf pzz.pdf temp.pdf press.pdf diff.pdf
DATA = $(shell ls data-*)
DIFF1 = $(subst data, diff1, $(DATA))
DIFF2 = $(subst data, diff2, $(DATA))
DIFF3 = $(subst data, diff3, $(DATA))

all: $(FIGURES) $(DIFF1) $(DIFF2) $(DIFF3)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: %.sh
	./$<

diff1-%.txt: asym-%.txt
	./compare.py $* data asym diff1

diff2-%.txt: dsmc-%.txt
	./compare.py $* dsmc data diff2

diff3-%.txt: dsmc-%.txt
	./compare.py $* dsmc asym diff3

clean:
	rm -f *{.tex,.aux,.log,-inc*} $(FIGURES)
