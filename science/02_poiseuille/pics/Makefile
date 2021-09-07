FIGURES = graph.pdf

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: %.sh
	./$<

clean:
	rm -f $(FIGURES)
