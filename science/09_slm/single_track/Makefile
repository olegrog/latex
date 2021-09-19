#VTK = $(shell ls *.vtk)
#PDF = $(VTK:.vtk=.pdf)

all: #$(PDF)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: _%.sh
	./$<

_%.sh: create_plots.py template.sh %.txt
	./$< $(word 2,$^) --name=$*

%.txt: extract_from_vtk.py %.vtk
	./$^ > $@

clean:
	rm -f $(PDF)
	rm -f *.{eps,log,tex,aux,txt}
	rm -f _*.sh

