PLOTS = width-speed width-specific_power depth-speed depth-specific_power
FIGURES = $(addsuffix .pdf,$(PLOTS))
SOURCE = $(wildcard *.txt)

options = \
    --xlabel=$(label_$2) \
    --ylabel=$(label_$1) \
    --base=$(formula_$2) \
    --func=$(formula_$1) \
    --c_num=$(column_numer_$1) \
    --c_exp=$(column_exper_$1) \
	--key=$(call key,$2)

label_power = 'P, W'
label_specific_power = 'P/U, \mathrm{J/mm}'
label_speed = 'U, mm/s'
label_width = 'd_\mathrm{track}-d_\mathrm{beam}, \mathrm{\mu m}'
label_depth = 'h_\mathrm{track}-h_\mathrm{powder}, \mathrm{\mu m}'

formula_speed = 'speed'
formula_specific_power = 'power/speed'
formula_width = 'value - 55'
formula_depth = 'value - 30'

column_exper_width = 3
column_numer_width = 2
column_exper_depth = 4
column_numer_depth = 3

key = $(if $(findstring specific,$1),'left','right')

all: $(FIGURES)

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: _%.sh $(SOURCE)
	./$<

_%.sh: create_plots.py template.sh
	./$^ --name=$* $(call options,$(word 1,$(subst -, ,$*)),$(word 2,$(subst -, ,$*)))

clean:
	rm -f $(FIGURES)
	rm -f *.{eps,log,tex,aux}
	rm -f _*.sh

