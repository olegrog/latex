UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	TIME = /usr/bin/time
endif
ifeq ($(UNAME), Darwin)
	TIME = gtime
endif
NUMBERS = $(shell python -c "import numpy, sys; numpy.savetxt(sys.stdout, [int(n) for n in numpy.logspace(.6,2.9,15)], fmt='%03d')")
MODELS = dvm lbm
FILES = $(filter-out %dvm-lbm.log, $(foreach m1,$(MODELS), $(addsuffix -$(m1).log, $(foreach m2,$(MODELS), $(addsuffix -$(m2), $(NUMBERS))))))
param = $(word $2, $(subst -, ,$1))

%.log:
	$(TIME) -v ../solver.py -e=50 -p=0 -m1=$(call param,$*,2) -m2=$(call param,$*,3) -N1=$(call param,$*,1)0 \
		-N2=$(call param,$*,1) -k=0.02 -r=5 > /dev/null 2> $@

all: $(FILES)
