# Copyright 2009--2010  Ed Bueler

all: lecture.pdf exers.pdf refs.pdf

zips: all
	rm -rf bueler_karthaus/  # clean out if prior existence
	mkdir bueler_karthaus/
	cp lecture.pdf exers.pdf refs.pdf bueler_karthaus/
	cp README bueler_karthaus/
	svn export mfiles/ bueler_karthaus/mfiles/
	(cd bueler_karthaus/mfiles/ && rm -rf nonpublic/)
	svn export petsc/ bueler_karthaus/petsc/
	(cd bueler_karthaus/petsc/ && cp notes/jacobiannotes.pdf .)
	(cd bueler_karthaus/petsc/ && rm -rf notes/)
	zip -r bueler_karthaus.zip bueler_karthaus/
	tar -cvzf bueler_karthaus.tar.gz bueler_karthaus/*
	rm -rf bueler_karthaus/

# list file names
figures := ice_shelf_edge_hg.jpg polarbear.jpg flowline.png fofv.png \
	palmer_land.png polaris.png schoof_planform.png schoof_sliders.png \
	streamisbrae.png siaerror.png Serac2.jpg eismintone.png \
	freehutter.png capnonflatobs.pdf classicalobs.pdf \
	brownian.pdf polythermal_types.pdf convanalysis.pdf \
	PISM_ross_speeds.png eisIIF.png earthcompare.png \
	NEgreenlandJoughin.png athabasca_cross.png hierarchy.pdf \
	joughin.png g3km_3_10_98.png g3km_3_10_98_hist.png \
	slabfigs.pdf  slabmasscontfig.pdf athabasca_deform.pdf \
	shelfthk.pdf shelfvel.pdf slabvel.pdf shelfconv.pdf shelfnumsoln.pdf \
	petscwww.pdf antinitial.pdf antfinal.pdf antvol.pdf green_transect.pdf \
	initialheat.pdf finalheat.pdf stability.pdf instability.pdf
epsfigures := diffstencil mahaffystencil sshape heatscaling siascaling \
        expstencil impstencil cnstencil exp2dstencil AofT slab
svgfigures := coffee heatconduction Rocket_nozzle_expansion
minputs := heat heatadapt halfar diffusion siaflat flowline \
	convanalysis ssaflowline ssainit testshelf verifysia
texinputs := intro.tex heatnum.tex sia.tex freebdry.tex ssa.tex next.tex

# this is just the first in the list of generated files for the animation:
animfigures := anim/heatmelt0.png anim/halfar0.png

# process file names
figures := $(addprefix photos/, $(figures))
epsfigures := $(addsuffix .pdf, $(addprefix pdffigs/, $(epsfigures)))
svgfigures := $(addsuffix .pdf, $(addprefix pdffigs/, $(svgfigures)))
minputs := $(addsuffix .m, $(addprefix mfiles/, $(minputs)))
mslim := $(subst .m,.slim.m,$(minputs))

# presentation
lecture.aux:  lecture.tex $(texinputs) \
	      $(epsfigures) $(svgfigures) $(figures) $(animfigures) $(mslim)
	pdflatex lecture
lecture.pdf:  lecture.aux
	pdflatex lecture

# exercises to accompany
exers.pdf:  exers.tex
	pdflatex exers
	pdflatex exers

# references to accompany; to create capture_bbl do this by hand:
#    make            # <-- creates lecture.aux
#    bibtex lecture
#    cp lecture.bbl capture_bbl  # <-- overwrites svn-controlled file
refs.pdf:  refs.tex capture_bbl
	pdflatex refs

pdffigs/%.pdf: pdffigs/%.eps
	epstopdf $< --outfile=$@

pdffigs/%.pdf: pdffigs/%.svg
	inkscape --export-eps=$*.eps $<
	epstopdf $*.eps --outfile=$@
	rm $*.eps

# see advice at http://www.ipgp.fr/~lucas/Contrib/animbeamer.html:
anim/heatmelt0.png: anim/Heat_eqn.gif
	(cd anim/ && convert Heat_eqn.gif heatmelt%d.png)
# change to 'matlab' if using that to generate movies
ANIMMATOCT := octave
anim/halfar0.png: mfiles/halfarmovie.m
	(cd mfiles/ && $(ANIMMATOCT) halfarmovie.m)
	(cd anim && mogrify -trim halfar*.png)

# the purpose of this nonsense is to put only the "meat" of codes into the
# PDF lectures.  thus we remove comments and print statements.  the flag
# "STRIPFROMHERE" can be used in a comment to remove everything after the flag
mfiles/%.slim.m: mfiles/%.m
	sed -e '/STRIPFROMHERE/,$$ d' $< > $@.tmp
	sed -e 's/ *fprintf/%/' -e 's/ *%/%/' -e '/^%/ d' -e 's/%.*//' $@.tmp > $@
	rm -rf $@.tmp

.PHONY: zips clean showEnv

showEnv:
	@echo 'figures = ' $(figures)
	@echo 'epsfigures = ' $(epsfigures)
	@echo 'svgfigures = ' $(svgfigures)
	@echo 'animfigures = ' $(animfigures)
	@echo 'minputs = ' $(minputs)
	@echo 'mslim = ' $(mslim)
	@echo 'texinputs = ' $(texinputs)

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.snm *.toc *.nav *.vrb \
	       mfiles/*.slim.m pdffigs/*.pdf anim/*.png bueler_karthaus.*

