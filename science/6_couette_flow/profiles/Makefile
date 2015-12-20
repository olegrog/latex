FIGURES = Pxy.pdf vx.pdf qx.pdf qy.pdf Pxx.pdf Pyy.pdf Pzz.pdf tau.pdf P.pdf omega.pdf

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: _%.sh
	./$<

_%.sh: template.sh plot_profiles.py
	./plot_profiles.py $*

clean:
	rm -f $(FIGURES)
