PLOTS = width-heat_transfer width-trumpf depth-heat_transfer
FIGURES = $(addsuffix .pdf,$(PLOTS))
SOURCE = $(wildcard *.txt)

options = \
    --column=$(column_$1) \
    --ground=$(ground_$1)

column_width = 3
column_depth = 4
ground_width = 55
ground_depth = 30

all: $(FIGURES)

.SECONDEXPANSION:

%.pdf: contour_plot.py ../$$(word 2,$$(subst -, ,$$*)).txt
	file=$(word 2,$^); \
	total=$$(awk '!/^#/{print $$1}' $$file | wc -l); \
	powers=$$(awk '!/^#/{print $$1}' $$file | uniq | wc -l); \
	./$^ $@ $$(($$total/$$powers)) $$powers -x=y -y=x --xlabel='U(mm/s)' --ylabel='P(W)' --logx --zoom=5 \
		$(call options,$(word 1,$(subst -, ,$*)))

../%.txt: ../%-*W.txt
	: > $@
	for f in $^; do s=$${f%%W.txt}; awk '!/^#/{ print '$${s##*-}', $$0 }' $$f >> $@; done

clean:
	rm -f $(FIGURES)

