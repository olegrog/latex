PROJECT=limiters

all: $(PROJECT).pdf

$(PROJECT).pdf: $(PROJECT).tex
	latexmk -pdf -shell-escape $(PROJECT)

clean:
	latexmk -C $(PROJECT)
