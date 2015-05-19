FIGURES := $(shell grep 'pdf:' Makefile | sed 's/:.*//' | grep -v FIGURES)
LOGS := $(shell grep '^log.p' Makefile | sed 's/:.*//')
DATADIR = ~/kesolver/data/couette
DIR = $(DATADIR)2.0
DIR2 = $(DATADIR)5.0

all: $(FIGURES) log.velocity_grid

kn0.1-center.pdf: $(DIR)/0.100/couette.kei $(DIR)/0.100/result/f24000.bin
	./out2f.py $? $@ center 10e-5 > /dev/null

kn0.1-boundary.pdf: $(DIR)/0.100/couette.kei $(DIR)/0.100/result/f24000.bin
	./out2f.py $? $@ boundary 10e-5 > /dev/null

kn1.0-center.pdf: $(DIR)/1.000/couette.kei $(DIR)/1.000/result/f6000.bin
	./out2f.py $? $@ center 30e-5 > /dev/null

kn1.0-boundary.pdf: $(DIR)/1.000/couette.kei $(DIR)/1.000/result/f6000.bin
	./out2f.py $? $@ boundary 10e-5 > /dev/null

kn10-center.pdf: $(DIR)/10.000/couette.kei $(DIR)/10.000/result/f18000.bin
	./out2f.py $? $@ center 30e-5 > /dev/null

kn10-boundary.pdf: $(DIR)/10.000/couette.kei $(DIR)/10.000/result/f18000.bin
	./out2f.py $? $@ boundary 10e-5 > /dev/null

kn1.0-center-xz5.pdf~: $(DIR2)/1.000/couette.kei $(DIR2)/1.000/result/f6000.bin
	./out2f.py $? $@ center 1.5e-5 y 9 8 > /dev/null

log.params-5.0:
	./discrete_error.py 2.5 20 26 4.3 1.27 > $@

log.params-2.0:
	./discrete_error.py 1.0 16 24 4.3 1.28 > $@

log.params-1.0:
	./discrete_error.py 0.5 14 23 4.3 1.29 > $@

log.params-0.1:
	./discrete_error.py .05 12 26 4.3 1.28 > $@

log.velocity_grid: $(LOGS)
	./collect_summary.sh > $@

clean:
	rm -f $(FIGURES)
