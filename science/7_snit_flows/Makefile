SHELL = /bin/bash
PROJECT = snit_flows
SUBDIRS = contours vs_kn vs_kn2 tikz
SOURCES = $(PROJECT).tex $(PROJECT).bib

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latex $(PROJECT)
	biber $(PROJECT)
	latex $(PROJECT)

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf}

m1:
	./discrete_error.py --temp-ratio=1.5 --radius=4.25 --velocity=1e-2 \
		--rN=8 --rtype=uniform --yN=8 --ytype=uniform

m2:
	./discrete_error.py --temp-ratio=1.5 --radius=5.3 --velocity=1e-2 \
		--rN=11 --rtype=hermite \
		--yN=26 --ymin=1.24 --ytype=geometric

m3:
	./discrete_error.py --temp-ratio=1.5 --radius=4.5 --velocity=1e-2 \
		--rN=12 --rmin=0.05  --rtype=quadratic \
		--yN=15 --ymin=0.005 --ytype=quadratic

m4:
	./discrete_error.py --temp-ratio=5 --radius=8 --velocity=1e-3 \
		--rN=16 --rmin=0.1 --rtype=quadratic \
		--yN=16 --ymin=0.1 --ytype=quadratic
	

