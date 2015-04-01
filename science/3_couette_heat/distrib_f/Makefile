FIGURES = kn0.1-center.pdf kn0.1-boundary.pdf kn1.0-center.pdf kn1.0-boundary.pdf kn10-center.pdf kn10-boundary.pdf
DATADIR = ~/kesolver/data/couette
U = 2.0
DIR = $(DATADIR)$(U)

all: $(FIGURES)

kn0.1-center.pdf: $(DIR)/0.100/couette.kei $(DIR)/0.100/result/f24000.bin
	./out2f.py $? 47 $@ 5e-5 > /dev/null

kn0.1-boundary.pdf: $(DIR)/0.100/couette.kei $(DIR)/0.100/result/f24000.bin
	./out2f.py $? 0 $@ 5e-5 > /dev/null

kn1.0-center.pdf: $(DIR)/1.000/couette.kei $(DIR)/1.000/result/f6000.bin
	./out2f.py $? 30 $@ 9e-5 > /dev/null

kn1.0-boundary.pdf: $(DIR)/1.000/couette.kei $(DIR)/1.000/result/f6000.bin
	./out2f.py $? 0 $@ 5e-5 > /dev/null

kn10-center.pdf: $(DIR)/10.000/couette.kei $(DIR)/10.000/result/f18000.bin
	./out2f.py $? 29 $@ 30e-5 > /dev/null

kn10-boundary.pdf: $(DIR)/10.000/couette.kei $(DIR)/10.000/result/f18000.bin
	./out2f.py $? 0 $@ 40e-5 > /dev/null

clean:
	rm -f $(FIGURES)
