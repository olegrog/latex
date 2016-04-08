MACRO = Pxy vx qx qy Pxx Pyy Pzz tau P omega 
FIGURES = $(addsuffix .pdf, $(MACRO))
CSV=$(addsuffix .csv, $(MACRO))

all: $(FIGURES) $(CSV)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: _%.sh
	./$<

_%.sh: template.sh plot_profiles.py
	./plot_profiles.py $*

%.csv: create_csv.py
	./$< $* | sed 's/+00/-00/g' > $@

clean:
	rm -f $(FIGURES)
