all: upwindnotes.pdf

codes   := 
figures := 

codes   := $(addprefix ../../mfiles/, $(codes))
figures := $(addprefix ../../photos/, $(figures))

# link ice-bib.bib needed for bibtex to work  (to file in pism-dev/doc/)

%.pdf: %.tex ice-bib.bib
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.toc *.snm *.bbl *.blg *.synctex.gz *.dvi *~

