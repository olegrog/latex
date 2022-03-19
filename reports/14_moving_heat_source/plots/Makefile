SCRIPT = pool.py
PLOTS = profiles
FIGURES = $(addsuffix .pdf,$(PLOTS))

all: $(FIGURES)

profiles.pdf: $(SCRIPT)
	./$< --Tratio=0.0378 --pdf -v --output=$@

.PHONY: all clean

clean:
	rm -f $(FIGURES)
