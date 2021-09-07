FIGURES = mean.pdf max.pdf rms.pdf

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: %.sh exclusions-*.txt
	sed -i~ 's/<file>/$*/' $<
	rm $<~
	chmod +x $< 
	./$<

mean.sh: template2.sh
	sed 's/<column>/2/;s/<operator>/mean/' $< > $@

max.sh: template2.sh
	sed 's/<column>/3/;s/<operator>/max/' $< > $@

rms.sh: template2.sh
	sed 's/<column>/4/;s/<operator>/rms/' $< > $@

clean:
	rm -f *{.tex,.aux,.log,-inc*} $(FIGURES)
