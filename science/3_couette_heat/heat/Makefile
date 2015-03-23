FIGURES = graph.pdf

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null

%.tex: %.sh
	./$<

clean:
	rm -f graph{.tex,.pdf,.aux,.log,-inc*}
