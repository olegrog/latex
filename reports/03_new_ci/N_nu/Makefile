PROJECT = $(notdir $(CURDIR))

all: $(PROJECT).pdf

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: %.sh *.txt
	./$<

.PHONY: clean

clean:
	rm -f $(PROJECT).pdf

