FILE=manuscript

all: $(FILE).pdf

$(FILE).pdf: $(FILE).tex
	latexmk -pdf $(FILE)

clean:
	latexmk -C
