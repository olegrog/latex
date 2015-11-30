FIGURES = delta_p.pdf
TXT = data.txt asym.txt snif.txt

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: plot.sh $(TXT)
	./$<

clean:
	rm -f *{.tex,.aux,.log,-inc*} $(FIGURES)
