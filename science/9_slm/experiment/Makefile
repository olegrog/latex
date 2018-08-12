PLOTS = width.pdf
SOURCE = single_track.txt awesumm.txt

all: $(PLOTS)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: %.sh $(SOURCE)
	./$<

clean:
	rm -f $(PDF)
	rm -f *.{eps,log,tex,aux}

