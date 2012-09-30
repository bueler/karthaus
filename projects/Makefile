all: upwindnotes.pdf linkmfiles

codes   := 
figures := 

codes   := $(addprefix ../mfiles/, $(codes))
figures := $(addprefix ../photos/, $(figures))

# link ice_bib.bib needed for bibtex to work  (to file in pism-dev/doc/)

%.pdf: %.tex ice_bib.bib
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*

linkmfiles:
	ln -s -f ../mfiles/diffusion.m .
	ln -s -f ../mfiles/exactshelf.m .
	ln -s -f ../mfiles/siageneral.m .

.PHONY: clean linkmfiles

clean:
	@rm -f *.out *.aux *.log *.toc *.snm *.bbl *.blg *.synctex.gz *.dvi *~

