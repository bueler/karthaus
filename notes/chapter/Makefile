all: kartbook-elb.pdf

%.pdf: %.tex
	pdflatex $*
	pdflatex $*

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.toc *.snm *.bbl *.blg *.synctex.gz *.dvi *.idx *~

