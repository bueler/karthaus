# Copyright 2009--2014  Ed Bueler

# build slides.pdf in ../ before this, otherwise it will not work

all: extraslides.pdf

# presentation
extraslides.aux:  extraslides.tex ps-freebdry.tex  ps-oversold.tex  ps-techskills.tex
	pdflatex extraslides
extraslides.pdf:  extraslides.aux
	pdflatex extraslides

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.snm *.toc *.nav *.vrb

