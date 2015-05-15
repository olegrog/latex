FIGURES = shear.pdf flow.pdf qflow.pdf qflowy.pdf pxx.pdf pzz.pdf temp.pdf press.pdf

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: %.sh
	./$<

clean:
	rm -f *{.tex,.aux,.log,-inc*} $(FIGURES)
