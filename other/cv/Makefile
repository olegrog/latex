FILE=cv

all: $(FILE).pdf

$(FILE).pdf: $(FILE).tex $(FILE).bib
	latex $(FILE)
	biber $(FILE)
	latex $(FILE)

