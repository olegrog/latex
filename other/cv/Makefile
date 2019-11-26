FILE = cv

all: $(FILE).pdf

$(FILE).pdf: $(FILE).tex main.bib other.bib photo.jpg
	latexmk -pdf $(FILE)

clean:
	rm $(FILE).{aux,bbl,bcf,blg,fdb_latexmk,fls,log,out,pdf,run.xml}
