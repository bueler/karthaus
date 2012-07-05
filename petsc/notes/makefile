all: jacobiannotes.pdf

# requires link for ice_bib.bib
jacobiannotes.pdf: jacobiannotes.tex
	pdflatex jacobiannotes.tex
	bibtex jacobiannotes
	pdflatex jacobiannotes.tex
	pdflatex jacobiannotes.tex

#### utilities

.PHONY : clean

clean :
	rm -f *.aux *.bbl *.blg *.log *.out

