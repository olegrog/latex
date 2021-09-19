all:
	latexmk -pdf -pdflatex='pdflatex --shell-escape'

clean:
	latexmk -C
	rm -f *.pgf-plot.{table,gnuplot} *.run.xml *.bbl

