PURE = pure.py
BINARY = binary.py

PDIR = pure
BDIR = binary

POPTS = -N=200 --pdf
BOPTS = -N=200 -K=0.5 --Dratio=0 --kratio=1 --pdf

PPLOTS = a0-V0.25 kV lambdaV
BPLOTS = bifurcation a0-G0.01-V0.2 lambda-G0.01-V0.2

PFIGURES = $(addprefix $(PDIR)/,$(addsuffix .pdf,$(PPLOTS)))
BFIGURES = $(addprefix $(BDIR)/,$(addsuffix .pdf,$(BPLOTS)))

parseGV = \
	-G=$(subst G,,$(word 1,$(subst -, ,$1))) \
	-V=$(subst V,,$(word 2,$(subst -, ,$1)))

all: mkdirs $(PFIGURES) $(BFIGURES)

$(PDIR)/a0-V%.pdf: $(PURE)
	./$< f $(POPTS) -V=$* -v -a --krange=0.05:2 --output=$@

$(PDIR)/kV.pdf: $(PURE)
	./$< 2 $(POPTS) -v -l -a --output=$@

$(PDIR)/lambdaV.pdf: $(PURE)
	./$< 2 $(POPTS) -v -l -a -w --output=$@

$(BDIR)/bifurcation.pdf: $(BINARY)
	./$< b $(BOPTS) --output=$@

$(BDIR)/a0-%.pdf: $(BINARY)
	./$< f $(BOPTS) $(call parseGV,$*) -v --krange=0.1:3 --output=$@

$(BDIR)/lambda-%.pdf: $(BINARY)
	./$< 2 $(BOPTS) -v -l -w --output=$@

.PHONY: all clean mkdirs

mkdirs:
	mkdir -p $(PDIR) $(BDIR)

clean:
	rm -f $(FIGURES)
