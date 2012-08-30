Build instructions
------

To build PDF form `slides.pdf` of the lectures, we require these tools

*  pdflatex
*  epstopdf
*  inkscape
*  pdftk
*  convert   (from ImageMagick)

and standard Linux tools including sed and GNU make.  See
Makefile for details.  The `pdffigs/` directory is filled by "make" 
below, from *.eps and *.svg sources, using epstopdf and inkscape.
Note that epstopdf is found in debian package texlive-extra-utils.

Octave and Matlab were used to generate some figures, typically via
.eps: `print -deps foo.eps`.

There is a ridiculously fragile route from .svg to get .pdf that will work
in pdflatex, using \includegraphics{}, without errors or warnings.
*But* at least it can be done at the command line.  Here is an example,
starting with an image "coffee.svg":

    $ inkscape --export-eps=coffee.eps coffee.svg
    $ epstopdf coffee.eps  # creates coffee.pdf which works in pdflatex
    $ rm coffee.eps

In this case coffee.pdf comes out *smaller* in file size than coffee.svg.
The disadvantage is that transparency is lost.  But at least there
are no errors in the .pdf that results from this route.

The Matlab codes in mfiles/ are stripped of comments and print
statements before going in the presentation.  That is the significance
of the *.slim.m names and the obscure sed command.

To clean up everything but leave the PDFs just produced:

    $ make clean

Building notes.pdf
-----

To build the notes in `notes/`, fewer tools are needed, but bibtex is required.
For BibTeX to work, a link must be created to the "master" .bib file from
the PISM source tree, or that file must be copied into the current directory.
Thus either do something similar to this:

    $ ln -s ~/pism-dev/doc/ice_bib.bib

Then do:

    $ make
