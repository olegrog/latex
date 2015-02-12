TRANSPORT=B.pdf B_4.pdf Q_2.pdf Q_3.pdf QQ_22.pdf QQ_3.pdf
NONLINEAR=J_1.pdf J_2.pdf J_3.pdf J_4.pdf J_5.pdf J_6.pdf
OTHER=diffB.pdf
SHELL=/bin/bash

all: $(TRANSPORT) $(NONLINEAR) $(OTHER)

B.pdf: ../tables/ci-functions.txt
	./plot.py $< 2 '\mathcal{B}' $@

B_4.pdf: ../data/gamma9.txt
	./plot.py $< 1 '\mathcal{B}^{(4)}' $@

Q_2.pdf: ../data/gamma8a.txt
	./plot.py $< 1 '\mathcal{Q}^{(0)}_2' $@

Q_3.pdf: ../data/gamma8a.txt
	./plot.py $< 2 '\mathcal{Q}^{(0)}_3' $@

QQ_22.pdf: ../data/gamma8b.txt
	./plot.py $< 2 '\tilde{\mathcal{Q}}^{(0)}_{22}' $@

QQ_3.pdf: ../data/gamma8b.txt
	./plot.py $< 3 '\tilde{\mathcal{Q}}^{(0)}_3' $@

J_1.pdf: ../data/gamma8ih.txt
	./plot.py $< 1 '\mathcal{J}^{(1)}_{x}' $@

J_2.pdf: ../data/gamma8ih.txt
	./plot.py $< 2 '\mathcal{J}^{(2)}_{xy}' $@

J_3.pdf: ../data/gamma8ih.txt
	./plot.py $< 3 '\mathcal{J}^{(3)}_{xy}' $@

J_4.pdf: ../data/gamma8ih.txt
	./plot.py $< 4 '\mathcal{J}^{(4)}_{xy}' $@

J_5.pdf: ../data/gamma8ih.txt
	./plot.py $< 5 '\mathcal{J}^{(5)}_{xyz}' $@

J_6.pdf: ../data/gamma8ih.txt
	./plot.py $< 6 '\mathcal{J}^{(6)}_{xyz}' $@

diffB.pdf: ../data/gamma1.txt
	./plot.py $< 3 '\displaystyle\frac{|\mathcal{B}-\mathcal{B}^*|}{\mathcal{B}^*}' $@ log .75

clean:
	rm -f $(TRANSPORT) $(NONLINEAR) $(OTHER)

