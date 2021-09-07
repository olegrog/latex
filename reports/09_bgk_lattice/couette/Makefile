TXT = $(shell find *.txt | grep -v ^k)
LOG = $(shell find log.*)
KNUDSENS = $(shell ls k*.txt | sed 's/.txt//;s/\./_/' | grep -oE '[0-9_]+' | sort | uniq)
FIGURES = $(TXT:.txt=.pdf) $(patsubst log.%,%.pdf,$(LOG)) $(patsubst %,shear-k%.pdf,$(KNUDSENS))
MODELS_FOR_SHEAR = $(shell grep suffix template-shear.sh | awk -F'[<"]' '{ print $$2 }')

options = \
	--buffer=$(call buffer,$1) \
	--kn=$(call knudsen,$1) \
	--qflow_factor=2 --shear_factor=0.5 \
	--hybrid=$(call hybrid,$1)

options_shear = \
	--buffer=$(call buffer,$1) \
	--kn=$(call knudsen,$1) \
	--exact=$(call exact,$1) \
	--suffix=$(call kn_suffix,$1) \
	--hybrid=1

options_solver = \
	--plot=0 $(call models,$1) $(call lattice,$1) $(call plot_norms,$1) \
	-N1=$(N1_$(call knudsen,$1)$(call extension,$1)) \
	-N2=$(N2_$(call knudsen,$1)$(call extension,$1)) \
	--refinement=$(refinement_$(call knudsen,$1)$(call extension,$1)) \
	--kn=$(call knudsen,$1) \
	--width=$(call width,$1) \
	--end=$(end_$(call knudsen,$1))


knudsen = $(or $(subst _,.,$(subst k,,$(filter k%,$1))),0.1)
extension = $(if $(findstring wide,$1),_wide,$(if $(findstring narrow,$1),_narrow))
width = $(width_$(call knudsen,$1)$(call extension,$1))
buffer = $(shell python -c "import math as m; print(0.5 - 2*$(call knudsen,$1)*$(call width,$1)/m.pi**.5)")
exact = $(shell tail -1 k$(call knudsen,$1).txt | awk '{ print $$3 }')
hybrid = $(if $(findstring hyb,$1),1,0)
kn_suffix = $(if $(filter-out 0.1,$(call knudsen,$1)),-$1)
models = $(if $(findstring dvm,$1), --model1=dvm,$(if $(findstring hyb,$1),, --model2=lbm))
lattice = $(addprefix  --lattice=,$(filter d3%,$1))
plot_norms = $(if $(findstring norms,$1), --plot-norms)
txt_for_shear = $(addsuffix $(call kn_suffix,$1).txt,$(MODELS_FOR_SHEAR))

N1_0.1 = 20
N1_0.1_wide = 8
N1_0.1_narrow = 22
N1_0.3 = 10
N1_0.03 = 20

N2_0.1 = 15
N2_0.1_wide = 20
N2_0.1_narrow = 8
N2_0.3 = 20
N2_0.03 = 12

refinement_0.1 = 5
refinement_0.1_wide = 5
refinement_0.1_narrow = 8
refinement_0.3 = 4
refinement_0.03 = 10

width_0.1 = 1.2
width_0.1_wide = 2.4
width_0.1_narrow = 0.6
width_0.3 = 0.8
width_0.03 = 3

end_0.1 = 20000
end_0.3 = 10000
end_0.03 = 30000

all: $(FIGURES) $(patsubst log.%,%.txt,$(LOG))

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: _%.sh
	./$<

_norms-%.sh: create_plots.py template-norms.sh norms-%.txt
	./$< $(word 2,$^) --name=norms-$* $(call options,$(subst -, ,$*))

_shear-%.sh: create_plots.py template-shear.sh $(call txt_for_shear,%)
	./$< $(word 2,$^) --name=shear-$* $(call options_shear,$*)

_%.sh: create_plots.py template.sh %.txt
	./$< $(word 2,$^) --name=$* $(call options,$(subst -, ,$*))

%.txt: log.%
	awk '/^#/{a=1} { if (a) print; }' $< > $@

log.%:
	./solver.py $(call options_solver,$(subst -, ,$*)) > $@

.PHONY: all clean $(SUBDIRS)

clean:
	rm -f $(FIGURES)
