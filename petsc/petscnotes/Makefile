all: jacobiannotes.pdf

# requires link for ice-bib.bib
jacobiannotes.pdf: jacobiannotes.tex
	pdflatex jacobiannotes.tex
	bibtex jacobiannotes
	pdflatex jacobiannotes.tex
	pdflatex jacobiannotes.tex

.PHONY : clean

clean :
	rm -f *.aux *.bbl *.blg *.log *.out

