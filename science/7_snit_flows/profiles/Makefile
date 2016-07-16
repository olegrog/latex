PATCHES = inner outer left bottom
SOLVERS = kes heat kgf asym-first asym-second asym-curv
FILES = $(foreach patch,$(PATCHES),$(addprefix $(patch)-, $(SOLVERS)))
FIGURES = $(addsuffix _temp.pdf, $(PATCHES)) $(addsuffix _vel.pdf, $(PATCHES))
VTK = ../elliptic/vtk
KN = 0.02
fields1 := T0,U1
fields2 := TMean,UMean
fields3 := T,U

all: $(FIGURES)

vtk_file = $(subst kgf,asym-first,$(VTK)/$1-$(KN).vtk)
fields = $(if $(filter heat kgf, $1),$(fields1),$(if $(filter kes, $1),$(fields2),$(fields3)))
for_kes = $(if $(filter kes, $1),--kn=$(KN))
corr = $(if $(filter heat kgf, $1),'x',$(if $(filter kes, $1),'x/np.sqrt(np.pi/2)/$(KN)','2*x/np.sqrt(np.pi)/$(KN)'))
options = $(call vtk_file,$1,$2) \
	--fields=$(call fields,$1) \
	--corr=$(call corr,$1) \
	--output=$2 \
	$(call for_kes,$1,$3)

xlabel = $(if $(filter left, $1),y,$(if $(filter bottom, $1),x,'\varphi'))
ylabel = $(if $(filter left bottom, $1),'T-T_\mathrm{HCE}','T-T_B')
xcoord = $(if $(filter left bottom, $1),0.31,0.28)

%.pdf: %.tex 
	pdflatex $< > /dev/null
	rm $*{-inc-eps-converted-to.pdf,-inc.eps,.aux,.log}

%.tex: _%.sh
	./$<

.SECONDEXPANSION:

inner-%.txt: patch_to_txt.py $$(call vtk_file,%)
	./$< $(call options,$*,$@) --equation='(x/0.3)**2 + (y/0.7)**2 - 1' --phi='np.arctan2(y,x)' --tau='[y*3/7,-x*7/3]' --exkind=quad

outer-%.txt: patch_to_txt.py $$(call vtk_file,%)
	./$< $(call options,$*,$@) --equation='(x/1.5)**2 + (y/1.0)**2 - 1' --phi='np.arctan2(y,x)' --tau='[-y*1.5,x/1.5]' --exkind=quad

left-%.txt: patch_to_txt.py $$(call vtk_file,%)
	./$< $(call options,$*,$@) --equation='x' --phi='y' --tau=[0,-1] $(if $(filter kes, $1),--exkind=xlog)

bottom-%.txt: patch_to_txt.py $$(call vtk_file,%)
	./$< $(call options,$*,$@) --equation='y' --phi='x' --tau=[1,0] $(if $(filter kes, $1),--exkind=xlog)

_%_temp.sh: create_plots.py template.sh $$(addprefix %-, $(addsuffix .txt, $(SOLVERS)))
	./$< $(word 2,$^) --name=$*_temp --patch=$* --column=2 --lmargin=4 \
		--ylabel=$(call ylabel,$*) --xlabel=$(call xlabel,$*) --xcoord=$(call xcoord,$*) --ytics=0.1 \
		--corr=1.0017

_%_vel.sh: create_plots.py template.sh $$(addprefix %-, $(addsuffix .txt, $(SOLVERS)))
	./$< $(word 2,$^) --name=$*_vel --patch=$* --column=3 --lmargin=5 \
		--ylabel='\frac{v_it_i}{k}' --xlabel=$(call xlabel,$*) --xcoord=0.28 --ytics=0.01

clean:
	rm -f $(FIGURES) rm *.{eps,log,tex}

.SECONDARY: $(addsuffix .txt, $(FILES))
