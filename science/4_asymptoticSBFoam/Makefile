FIGURES=transport/Y1_2.pdf tikz/geometry.pdf simple/residual.pdf \
	plane/T-heat.pdf plane/T-asym.pdf plane/U.pdf plane/U01.pdf \
	noncoaxial/cylinder5.pdf noncoaxial/sphere5.pdf noncoaxial/inverse5.pdf noncoaxial/cylinder2.pdf \
	force/forces.pdf force/forces-close.pdf force/profiles.pdf \
	force/terms-cylinder-inner.pdf force/terms-cylinder-inner-inv.pdf force/terms-cylinder-outer.pdf \
	elliptic/U.pdf elliptic/delta-x.pdf \
	moment/profiles.pdf moment/moment-beta.pdf \
	temper/F-tau-cylinders-inner.pdf temper/F-tau-cylinders-outer.pdf temper/U-tau-cylinders.pdf
SHELL=/bin/bash

all: $(FIGURES)

transport/Y1_2.pdf: transport/Y1_2.txt transport/Y1.py
	cd transport && ./Y1.py

tikz/geometry.pdf: tikz/geometry.tex
	cd tikz && latex geometry > log.latex

plane/T-heat.pdf: plane/heat.raw plane/temperature.py
	cd plane && ./temperature.py heat

plane/T-asym.pdf: plane/asym.raw plane/temperature.py
	cd plane && ./temperature.py asym

plane/U.pdf: plane/asym.raw plane/velocity.py
	cd plane && ./velocity.py 0

plane/U01.pdf: plane/asym.raw plane/velocity.py
	cd plane && ./velocity.py 0.01

noncoaxial/cylinder5.pdf: noncoaxial/cylinder5.raw noncoaxial/velocity.py
	cd noncoaxial && ./velocity.py cylinder5

noncoaxial/sphere5.pdf: noncoaxial/sphere5.raw noncoaxial/velocity.py
	cd noncoaxial && ./velocity.py sphere5

noncoaxial/inverse5.pdf: noncoaxial/inverse5.raw noncoaxial/velocity.py
	cd noncoaxial && ./velocity.py inverse5

noncoaxial/cylinder2.pdf: noncoaxial/cylinder2.raw noncoaxial/velocity.py
	cd noncoaxial && ./velocity.py cylinder2

force/forces.pdf: force/spheres.txt force/cylinders.txt force/force.py
	cd force && ./force.py

force/forces-close.pdf: force/spheres.txt force/cylinders.txt force/force-close.py
	cd force && ./force-close.py

force/profiles.pdf: force/*Force.raw force/profile.py
	cd force && ./profile.py

force/terms-cylinder-inner.pdf: force/cylinder-inner[HVT]*.raw force/terms.py
	cd force && ./terms.py cylinder-inner 0.5

force/terms-cylinder-inner-inv.pdf: force/cylinder-inner-inv[HVT]*.raw force/terms.py
	cd force && ./terms.py cylinder-inner-inv 0.5

force/terms-cylinder-outer.pdf: force/cylinder-outer[HVT]*.raw force/terms.py
	cd force && ./terms.py cylinder-outer 0

elliptic/U.pdf: elliptic/elliptic.raw elliptic/rotated.raw elliptic/velocity.py
	cd elliptic && ./velocity.py

elliptic/delta-x.pdf: elliptic/heat5.txt elliptic/asym5.txt elliptic/delta.py
	cd elliptic && ./delta.py

moment/profiles.pdf: moment/beta*.raw moment/profile.py
	cd moment && ./profile.py

moment/moment-beta.pdf: moment/elliptic.txt moment/inverse.txt moment/moment.py
	cd moment && ./moment.py

temper/U-tau-cylinders.pdf: temper/cylinders.txt temper/U-tau.py
	cd temper && ./U-tau.py cylinders

temper/F-tau-cylinders-inner.pdf: temper/cylinders.txt temper/F-tau.py
	cd temper && ./F-tau.py cylinders inner

temper/F-tau-cylinders-outer.pdf: temper/cylinders.txt temper/F-tau.py
	cd temper && ./F-tau.py cylinders outer

simple/residual.pdf: simple/convergence.py simple/p_0  simple/T_0  simple/Ux_0  simple/Uy_0
	cd simple && ./convergence.py

clean:
	rm -f $(FIGURES)

