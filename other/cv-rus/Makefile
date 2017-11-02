FILE = cv-rus

all: $(FILE).pdf

$(FILE).pdf: $(FILE).tex main.bib other.bib
	latexmk -pdf $(FILE)

clean:
	rm $(FILE).{aux,bbl,bcf,blg,fdb_latexmk,fls,log,out,pdf,run.xml}
