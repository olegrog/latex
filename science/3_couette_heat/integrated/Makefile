FIGURES = shear.pdf flow.pdf qflow.pdf qflowy.pdf pxx.pdf pzz.pdf temp.pdf

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null

%.tex: %.sh
	./$<

clean:
	rm -f {shear,flow,qflow,qflowy,pxx,pzz,temp}{.tex,.pdf,.aux,.log,-inc*}
