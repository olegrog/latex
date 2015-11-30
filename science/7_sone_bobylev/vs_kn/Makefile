FIGURES = top_temp.pdf top_flow.pdf bottom_temp.pdf bottom_flow.pdf
SCRIPTS=$(FIGURES:.pdf=.sh)
TXT = data.txt asym.txt snif.txt heat.txt high.txt

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: %.sh $(TXT)
	./$<

$(SCRIPTS): create_plots.py template.sh
	./$< $@

clean:
	rm -f *{.tex,.aux,.log,-inc*} $(SCRIPTS) $(FIGURES)
