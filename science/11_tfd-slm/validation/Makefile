PLOTS = width-specific_power depth-specific_power
FIGURES = $(addsuffix .pdf,$(PLOTS))
SOURCE = $(wildcard *.txt)

options = \
    --xlabel=$(label_$2) \
    --ylabel=$(label_$1) \
    -x=$(name_$2) \
    -y=$1 \
    --logx --use-tex

label_power = '$$P$$, W'
label_specific_power = '$$P/U$$, J/mm'
label_speed = '$$U$$, mm/s'
label_width = '$$d$$, µm'
label_depth = '$$h$$, µm'

name_speed = 'speed'
name_power = 'power'
name_specific_power = 'power/speed'

all: $(FIGURES)

%.pdf: plot.py $(SOURCE)
	./$< $@ $(call options,$(word 1,$(subst -, ,$*)),$(word 2,$(subst -, ,$*)))

clean:
	rm -f $(FIGURES)

