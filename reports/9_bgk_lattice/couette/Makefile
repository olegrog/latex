MODELS = $(shell ls *.txt | grep -v ^k)
FIGURES = $(MODELS:.txt=.pdf)

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: _%.sh
	./$<

_%.sh: create_plots.py template.sh %.txt
	./$< $(word 2,$^) --name=$*

clean:
	rm -f $(FIGURES)
