MODELS = $(shell ls *.txt | grep -v ^k)
FIGURES = $(MODELS:.txt=.pdf)

buffer = $(if $(filter %wide, $1),'0.25',$(if $(filter %narrow, $1),'0.43229725','0.3645945'))

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: _%.sh
	./$<

_norms-%.sh: create_plots.py template-norms.sh norms-%.txt
	./$< $(word 2,$^) --name=norms-$* --buffer=$(call buffer,$*)

_%.sh: create_plots.py template.sh %.txt
	./$< $(word 2,$^) --name=$* --buffer=$(call buffer,$*)

clean:
	rm -f $(FIGURES)
