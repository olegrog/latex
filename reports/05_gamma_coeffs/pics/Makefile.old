FIRST=A.pdf
SECOND=B.pdf B_4.pdf
THIRD=T1_1.pdf T1_2.pdf T2_1.pdf T2_2.pdf TT_12.pdf TT_2.pdf
FOURTH=Q_2.pdf Q_3.pdf QQ_22.pdf QQ_3.pdf
TRANSPORT=$(FIRST) $(SECOND) $(THIRD) $(FOURTH)
NONLINEAR=J_a.pdf J_b.pdf J_1.pdf J_2.pdf J_3.pdf J_4.pdf J_5.pdf J_6.pdf
OTHER=diffB.pdf
SHELL=/bin/bash

all: $(TRANSPORT) $(NONLINEAR) $(OTHER)

A.pdf: ../tables/ci-functions.txt
	./plot.py $< 1 '\mathcal{A}' $@

B.pdf: ../tables/ci-functions.txt
	./plot.py $< 2 '\mathcal{B}' $@

B_4.pdf: ../data/gamma9.txt
	./plot.py $< 1 '\mathcal{B}^{(4)}' $@

T1_1.pdf: ../data/gamma10a.txt
	./plot.py $< 1 '\mathcal{T}^{(1)}_1' $@

T1_2.pdf: ../data/gamma10a.txt
	./plot.py $< 2 '\mathcal{T}^{(1)}_2' $@

T2_1.pdf: ../data/gamma10b.txt
	./plot.py $< 1 '\mathcal{T}^{(2)}_1' $@

T2_2.pdf: ../data/gamma10b.txt
	./plot.py $< 2 '\mathcal{T}^{(2)}_2' $@

TT_12.pdf: ../data/gamma10c.txt
	./plot.py $< 1 '\tilde{\mathcal{T}}^{(0)}_{12}' $@ none default 2.5

TT_2.pdf: ../data/gamma10c.txt
	./plot.py $< 2 '\tilde{\mathcal{T}}^{(0)}_2' $@

Q_2.pdf: ../data/gamma8a.txt
	./plot.py $< 1 '\mathcal{Q}^{(0)}_2' $@

Q_3.pdf: ../data/gamma8a.txt
	./plot.py $< 2 '\mathcal{Q}^{(0)}_3' $@

QQ_22.pdf: ../data/gamma8b.txt
	./plot.py $< 2 '\tilde{\mathcal{Q}}^{(0)}_{22}' $@

QQ_3.pdf: ../data/gamma8b.txt
	./plot.py $< 3 '\tilde{\mathcal{Q}}^{(0)}_3' $@

J_a.pdf: ../data/gamma10ih.txt
	./plot.py $< 1 '\mathcal{J}^{(\mathrm{I})}_{xy}' $@ none default 2.5

J_b.pdf: ../data/gamma10ih.txt
	./plot.py $< 2 '\mathcal{J}^{(\mathrm{II})}_{xyz}' $@ none default 2.5

J_1.pdf: ../data/gamma8ih.txt
	./plot.py $< 1 '\mathcal{J}^{(1)}_{x}' $@

J_2.pdf: ../data/gamma8ih.txt
	./plot.py $< 2 '\mathcal{J}^{(2)}_{xy}' $@ none default 2.5

J_3.pdf: ../data/gamma8ih.txt
	./plot.py $< 3 '\mathcal{J}^{(3)}_{xy}' $@

J_4.pdf: ../data/gamma8ih.txt
	./plot.py $< 4 '\mathcal{J}^{(4)}_{xy}' $@

J_5.pdf: ../data/gamma8ih.txt
	./plot.py $< 5 '\mathcal{J}^{(5)}_{xyz}' $@ none default 2.5

J_6.pdf: ../data/gamma8ih.txt
	./plot.py $< 6 '\mathcal{J}^{(6)}_{xyz}' $@

diffB.pdf: ../data/gamma1.txt
	./plot.py $< 3 '\displaystyle\frac{|\mathcal{B}-\mathcal{B}^*|}{\mathcal{B}^*}' $@ log .75

clean:
	rm -f $(TRANSPORT) $(NONLINEAR) $(OTHER)

