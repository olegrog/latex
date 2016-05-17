FIGURES = top_temp.pdf top_flow.pdf bottom_temp.pdf bottom_flow.pdf
TXT = data.txt high.txt asym1.txt asym2.txt

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: _%.sh $(TXT)
	./$<

_%.sh: create_plots.py template.sh
	./$< $*

clean:
	rm -f $(FIGURES)
