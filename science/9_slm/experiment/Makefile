PLOTS = width-speed width-specific_power
FIGURES = $(addsuffix .pdf,$(PLOTS))
SOURCE = $(wildcard *.txt)

options = \
    --xlabel=$(call label,$2) \
    --ylabel=$(call label,$1) \
    --base=$(call base,$2) \
	--key=$(call key,$2)

label = $(label_$1)
base = $(base_$1)
key = $(if $(findstring specific,$1),'left','right')

label_power = 'P, W'
label_specific_power = 'P/U, \mathrm{J/mm}'
label_speed = 'U, mm/s'
label_width = 'd_\mathrm{track}-d_\mathrm{beam}, \mathrm{\mu m}'

base_speed = speed
base_specific_power = power/speed

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

