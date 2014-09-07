FIGURES=transport/Y1_2.pdf tikz/geometry.pdf plane/T-heat.pdf plane/T-asym.pdf plane/U.pdf \
	plane/U01.pdf noncoaxial/cylinder5.pdf noncoaxial/sphere5.pdf noncoaxial/inverse5.pdf  \
	noncoaxial/cylinder2.pdf force/forces.pdf force/forces-inv.pdf force/forces-alpha.pdf  \
	force/profiles.pdf force/terms-cylinder.pdf force/terms-cylinder-inv.pdf force/terms-sphere.pdf \
	force/terms-sphere-inv.pdf elliptic/U.pdf elliptic/delta-x.pdf elliptic/maxU-alpha.pdf \
	moment/profiles.pdf moment/moment-beta.pdf simple/residual.pdf
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

force/forces-inv.pdf: force/spheres-inv.txt force/cylinders.txt force/force-inv.py
	cd force && ./force-inv.py

force/forces-alpha.pdf: force/alpha.txt force/force-alpha.py
	cd force && ./force-alpha.py

force/profiles.pdf: force/*Force.raw force/profile.py
	cd force && ./profile.py

force/terms-cylinder.pdf: force/cylinder[PS]*.raw force/terms-new.py
	cd force && ./terms-new.py cylinder 0.5

force/terms-cylinder-inv.pdf: force/cylinder-inv[PS]*.raw force/terms-new.py
	cd force && ./terms-new.py cylinder-inv 0.5

force/terms-sphere.pdf: force/sphere[PS]*.raw force/terms-new.py
	cd force && ./terms-new.py sphere 0.5

force/terms-sphere-inv.pdf: force/sphere-inv[PS]*.raw force/terms-new.py
	cd force && ./terms-new.py sphere-inv 0.5

elliptic/U.pdf: elliptic/elliptic.raw elliptic/rotated.raw elliptic/velocity.py
	cd elliptic && ./velocity.py

elliptic/delta-x.pdf: elliptic/heat5.txt elliptic/asym5.txt elliptic/delta.py
	cd elliptic && ./delta.py

elliptic/maxU-alpha.pdf: elliptic/alpha.txt elliptic/alpha.py
	cd elliptic && ./alpha.py

moment/profiles.pdf: moment/beta*.raw moment/profile.py
	cd moment && ./profile.py

moment/moment-beta.pdf: moment/elliptic.txt moment/inverse.txt moment/moment.py
	cd moment && ./moment.py

simple/residual.pdf: simple/convergence.py simple/p_0  simple/T_0  simple/Ux_0  simple/Uy_0
	cd simple && ./convergence.py

clean-pics:
	rm -f $(FIGURES)

