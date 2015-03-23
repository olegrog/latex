FIGURES = error.pdf

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null

%.tex: %.sh
	./$<

clean:
	rm -f error{.tex,.pdf,.aux,.log,-inc*}
