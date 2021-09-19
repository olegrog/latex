PLOTS = convergence korobov
FIGURES = $(addsuffix .pdf,$(PLOTS))

all: $(FIGURES)

convergence.pdf: plot.py *.txt
	./$<

korobov.pdf: plot2.py *.txt
	./$<

.PHONY: clean

clean:
	rm -f $(FIGURES)

