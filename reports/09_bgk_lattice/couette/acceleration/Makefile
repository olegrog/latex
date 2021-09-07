FIGURES = acceleration.pdf

all: $(FIGURES)

%.pdf: %.tex *.txt
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: _%.sh
	./$<

_%.sh: %.sh
	sed 's/<name>/'$*'/' $^ > $@
	chmod +x $@

clean:
	rm -f $(FIGURES)
