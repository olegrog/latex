PDF = width.pdf
SOURCE = $(shell ls *.txt)

all: $(PDF)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: %.sh $(SOURCE)
	./$<

clean:
	rm -f $(PDF)
	rm -f *.{eps,log,tex,aux}

