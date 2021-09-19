PLOT = convergence

all: $(PLOT).pdf

%.pdf: plot.py *.txt
	./$<

.PHONY: clean

clean:
	rm -f $(PLOT).pdf

