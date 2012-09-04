all: lecture.pdf

# need to do "make" in ../ for this too work
mc_zip: lecture.pdf
	rm -rf mc_zip/
	svn export ../mfiles/ mc_zip/
	cp lecture.pdf mc_zip/
	cp refs_stable.pdf mc_zip/refs.pdf

mccarthy-bueler.zip: mc_zip
	zip -r mccarthy-bueler.zip mc_zip/

mccarthy-bueler.tar.gz: mc_zip
	tar -cvzf mccarthy-bueler.tar.gz mc_zip/*

# presentation; need to do "make" in ../ for this too work
lecture.pdf:  lecture.tex style.tex preamble.tex sia.tex masscont.tex ssa.tex advanced.tex
	pdflatex lecture
	pdflatex lecture

# probably need to uncomment lines in lecture.tex to generate refs
refs_stable.pdf: lecture.pdf refs.tex
	bibtex lecture
	pdflatex lecture
	pdflatex lecture
	cp lecture.bbl capture_bbl.tex
	pdflatex refs
	pdflatex refs
	cp refs.pdf refs_stable.pdf

# a good way to get a copy for print with multiple pages per sheet is:
printable.pdf:
	pdfnup --nup '1x2' --scale 0.9 --frame true --no-landscape -- lecture.pdf -
	pdftk lecture-nup.pdf refs_stable.pdf cat output printable.pdf

.PHONY: clean mc_zip

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.snm *.toc *.nav *.vrb *.synctex.gz
	@rm -rf mc_zip/

