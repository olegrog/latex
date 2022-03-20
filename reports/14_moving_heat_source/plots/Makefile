SCRIPT = pool.py
PLOTS = boundary gradient
OPTS = --Tratio=0.0378 --pdf -v -N=200
FIGURES = $(addsuffix .pdf,$(PLOTS))

all: $(FIGURES)

boundary.pdf: $(SCRIPT)
	./$< b $(OPTS) --output=$@

gradient.pdf: $(SCRIPT)
	./$< g $(OPTS) --output=$@ --figsize=5:3

.PHONY: all clean

clean:
	rm -f $(FIGURES)
