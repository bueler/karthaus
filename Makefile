# Copyright 2009--2012  Ed Bueler

all: slides.pdf

# list file names
figures := iceshelfedge.jpg flowline.png fofv.png \
	palmer_land.png polaris.png schoof_planform.png schoof_sliders.png \
	streamisbrae.png siaerror.png Serac2.jpg eismintone.png \
	freehutter.png capnonflatobs.pdf classicalobs.pdf brownian.pdf \
	polythermal_types.pdf convanalysis.pdf eisIIF.png earthcompare.png \
	NEgreenlandJoughin.png athabasca_cross.png hierarchy.pdf \
	joughin.png g3km_3_10_98.png g3km_3_10_98_hist.png green_transect.pdf \
	slabfigs.pdf  slabmasscontfig.pdf athabasca_deform.pdf \
	shelfthk.pdf shelfvel.pdf slabvel.pdf shelfconv.pdf shelfnumsoln.pdf \
	antinitial.png antfinal.png antvolcompare.pdf antthickchange.png \
	initialheat.pdf finalheat.pdf stability.pdf instability.pdf slab.pdf \
	roughfinal.png roughinitial.png roughtimesteps.pdf ssavel8km.pdf \
	expstencil.pdf exp2dstencil.pdf diffstencil.pdf mahaffystencil.pdf \
	heatscaling.pdf
epsfigures :=  sshape siascaling impstencil cnstencil AofT
svgfigures := coffee heatconduction Rocket_nozzle_expansion
minputs := heat heatadapt halfar diffusion siaflat flowline \
	convanalysis ssaflowline ssainit testshelf verifysia
texinputs := intro.tex sia.tex masscont.tex ssa.tex

# this is just the first in the list of generated files for the animation:
animfigures := anim/heatmelt0.png anim/halfar0.png

# process file names
figures := $(addprefix photos/, $(figures))
epsfigures := $(addsuffix .pdf, $(addprefix pdffigs/, $(epsfigures)))
svgfigures := $(addsuffix .pdf, $(addprefix pdffigs/, $(svgfigures)))
minputs := $(addsuffix .m, $(addprefix mfiles/, $(minputs)))
mslim := $(subst .m,.slim.m,$(minputs))

# presentation
slides.aux:  slides.tex $(texinputs) \
	      $(epsfigures) $(svgfigures) $(figures) $(animfigures) $(mslim)
	pdflatex slides
slides.pdf:  slides.aux
	pdflatex slides

pdffigs/%.pdf: pdffigs/%.eps
	epstopdf $< --outfile=$@

pdffigs/%.pdf: pdffigs/%.svg
	inkscape --export-eps=$*.eps $<
	epstopdf $*.eps --outfile=$@
	rm $*.eps

# see advice at http://www.ipgp.fr/~lucas/Contrib/animbeamer.html:
anim/heatmelt0.png: anim/Heat_eqn.gif
	(cd anim/ && convert Heat_eqn.gif heatmelt%d.png)

# the old way to generate the Halfar movie:
## change to 'matlab' if using that to generate movies
#ANIMMATOCT := octave
#anim/halfar0.png: mfiles/halfarmovie.m
#	(cd mfiles/ && $(ANIMMATOCT) halfarmovie.m)
#	(cd anim && mogrify -trim halfar*.png)

# the new way to generate the Halfar movie:
anim/halfar0.png: anim/halfarfigs.tar.gz
	(cd anim/ && tar -xzvf halfarfigs.tar.gz)

# the purpose of this nonsense is to put only the "meat" of codes into the
# PDF lectures.  thus we remove comments and print statements.  the flag
# "STRIPFROMHERE" can be used in a comment to remove everything after the flag
mfiles/%.slim.m: mfiles/%.m
	sed -e '/STRIPFROMHERE/,$$ d' $< > $@.tmp
	sed -e 's/ *fprintf/%/' -e 's/ *%/%/' -e '/^%/ d' -e 's/%.*//' $@.tmp > $@
	rm -rf $@.tmp

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.snm *.toc *.nav *.vrb \
	       mfiles/*.slim.m pdffigs/*.pdf anim/*.png

