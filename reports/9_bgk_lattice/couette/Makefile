SOURCES = $(shell ls *.txt | grep -v ^k)
KNUDSENS = $(shell ls k*.txt | sed 's/.txt//;s/\./_/' | grep -oE '[0-9_]+' | sort | uniq)
FIGURES = $(SOURCES:.txt=.pdf) $(patsubst %,shear-k%.pdf,$(KNUDSENS))

options = \
	--buffer=$(call buffer,$1) \
	--kn=$(call knudsen,$1) \
	--qflow=$(qflow_$(call knudsen,$1)) \
	--ymin=$(or $(ymin_$(call knudsen,$1)),0) \
	--hybrid=$(call hybrid,$1)

options_shear = \
	--buffer=$(call buffer,$1) \
	--exact=$(exact_$(call knudsen,$1)) \
	--suffix=$(call kn_suffix,$1) \
	--hybrid=1


extension = $(if $(findstring wide,$1),_wide,$(if $(findstring narrow,$1),_narrow))
buffer = $(buffer_$(call knudsen,$1)$(call extension,$1))
knudsen = $(or $(subst _,.,$(subst k,,$(filter k%,$1))),0.1)
hybrid = $(if $(findstring hyb,$1),1,0)
kn_suffix = $(if $(filter-out 0.1,$(call knudsen,$1)),-$1)

buffer_0.1 = 0.3646
buffer_0.1_wide = 0.25
buffer_0.1_narrow = 0.4323
buffer_0.3 = 0.25
buffer_0.03 = 0.398

qflow_0.1 = 40
qflow_0.3 = 10
qflow_0.03 = 80

exact_0.1 = -0.0415561
exact_0.3 = -0.0934498
exact_0.03 = -0.0141380

ymin_0.03 = -0.02

all: $(FIGURES)

subdirs:
	for dir in $(SUBDIRS); do \
		cp *.py Makefile $$dir; \
		$(MAKE) -C $$dir; \
		rm $$dir/{*.py,Makefile}; \
	done

%.pdf: %.tex
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: _%.sh
	./$<

_norms-%.sh: create_plots.py template-norms.sh norms-%.txt
	./$< $(word 2,$^) --name=norms-$* $(call options,$(subst -, ,$*))

_shear-%.sh: create_plots.py template-shear.sh
	./$< $(word 2,$^) --name=shear-$* $(call options_shear,$*)

_%.sh: create_plots.py template.sh %.txt
	./$< $(word 2,$^) --name=$* $(call options,$(subst -, ,$*))

clean:
	rm -f $(FIGURES)
