FIGURES = Pxy.pdf vx.pdf qx.pdf qy.pdf Pxx.pdf Pyy.pdf Pzz.pdf tau.pdf P.pdf omega.pdf
SCRIPTS=$(FIGURES:.pdf=.sh)

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: %.sh
	./$<

$(SCRIPTS): template.sh plot_profiles.py
	./plot_profiles.py $@

clean:
	rm -f *{.tex,.aux,.log,-inc*} $(SCRIPTS) $(FIGURES)
