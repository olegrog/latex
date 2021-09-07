SOURCES := $(wildcard *.tex)
OBJECTS := $(patsubst %.tex,%.pdf,$(SOURCES))

all: $(OBJECTS)

%.pdf: %.tex
	latex $* > log.$*

