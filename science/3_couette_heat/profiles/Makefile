FIGURES = Pxy.pdf vx.pdf qx.pdf qy.pdf Pxx.pdf Pzz.pdf tau.pdf
SCRIPTS=$(FIGURES:.pdf=.sh)

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*-inc-eps-converted-to.pdf

%.tex: %.sh
	./$<

$(SCRIPTS):
	./plot_profiles.py $@

clean:
	rm -f {Pxy,vx,qx,qy,Pxx,Pzz,tau}{.tex,.pdf,.aux,.log,-inc*,.sh}
