FILE=limiters

all: $(FILE).pdf

$(FILE).pdf: $(FILE).tex 
	latex $(FILE)
	latex --shell-escape $(FILE)

