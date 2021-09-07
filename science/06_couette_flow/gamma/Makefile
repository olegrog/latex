SECOND=B.pdf B_4.pdf
THIRD=T1_1.pdf T1_2.pdf T2_1.pdf T2_2.pdf TT_12.pdf TT_2.pdf
FOURTH=Q_2.pdf Q_3.pdf QQ_22.pdf QQ_3.pdf
TRANSPORT=$(SECOND) $(THIRD) $(FOURTH)
CSV=gamma1.csv gamma2.csv
DIR=../../../reports/5_gamma_coeffs/infty
SHELL=/bin/bash

all: $(TRANSPORT) $(CSV)

B.pdf: $(DIR)/B.txt plot.py
	./plot.py $< 1 '\mathcal{B}' $@

B_4.pdf: $(DIR)/B4.txt plot.py
	./plot.py $< 1 '\mathcal{B}^{(4)}' $@

T1_1.pdf: $(DIR)/T1_1.txt plot.py
	./plot.py $< 1 '\mathcal{T}^{(1)}_1' $@

T1_2.pdf: $(DIR)/T1_2.txt plot.py
	./plot.py $< 1 '\mathcal{T}^{(1)}_2' $@

T2_1.pdf: $(DIR)/T2_1.txt plot.py
	./plot.py $< 1 '\mathcal{T}^{(2)}_1' $@

T2_2.pdf: $(DIR)/T2_2.txt plot.py
	./plot.py $< 1 '\mathcal{T}^{(2)}_2' $@

TT_12.pdf: $(DIR)/TT12.txt plot.py
	./plot.py $< 1 '\tilde{\mathcal{T}}^{(0)}_{12}' $@ none default 2.5

TT_2.pdf: $(DIR)/TT2.txt plot.py
	./plot.py $< 1 '\tilde{\mathcal{T}}^{(0)}_2' $@

Q_2.pdf: $(DIR)/Q2.txt plot.py
	./plot.py $< 1 '\mathcal{Q}^{(0)}_2' $@

Q_3.pdf: $(DIR)/Q3.txt plot.py
	./plot.py $< 1 '\mathcal{Q}^{(0)}_3' $@

QQ_22.pdf: $(DIR)/QQ22.txt plot.py
	./plot.py $< 1 '\tilde{\mathcal{Q}}^{(0)}_{22}' $@

QQ_3.pdf: $(DIR)/QQ3.txt plot.py
	./plot.py $< 1 '\tilde{\mathcal{Q}}^{(0)}_3' $@

gamma1.csv: create_csv.py
	./$< $(DIR) 0 > $@

gamma2.csv: create_csv.py
	./$< $(DIR) 1 > $@

clean:
	rm -f $(TRANSPORT) $(NONLINEAR) $(OTHER)

