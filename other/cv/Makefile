FILE=cv

all: $(FILE).pdf

$(FILE).pdf: $(FILE).tex $(FILE).bib
	latex $(FILE)
	biber $(FILE)
	latex $(FILE)

clean:
	rm $(FILE).{aux,bbl,bcf,blg,fdb_latexmk,fls,log,out,pdf,run.xml}
