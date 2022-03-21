SCRIPT = pool.py
PLOTS = boundary gradient speed rate
OPTS = --Tratio=0.0378 --pdf -v -a -N=200 --figsize=5:3
FIGURES = $(addsuffix .pdf,$(PLOTS))

all: $(FIGURES)

%.pdf: $(SCRIPT)
	./$< $* $(OPTS) --output=$@

.PHONY: all clean

clean:
	rm -f $(FIGURES)
