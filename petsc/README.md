SSA solver using PETSc
========

Copyright 2009--2014  Ed Bueler

The SSA code `ssaflowline.c` solves the same problem solved by the Matlab/Octave
codes `testshelf.m`, `ssaflowline.m`, and `flowline.m` in directory `../mfiles/`.

See the PDF lecture `../slides.pdf` for a description of the problem and of the
Picard iteration solution implemented in the Matlab/Octave codes.

See `petscnotes/jacobiannotes.pdf` for the purpose of this example, and many
details of its construction.

To build and run in this directory:
  0)  needs PETSc version 3.5
  1)  make sure `PETSC_DIR` and `PETSC_ARCH` are set correctly
  2)  build:
        $ make
  3)  run with defaults:
        $ ./ssaflowline
  4)  get help and suggestions on other options combinations:
        $ ./ssaflowline -help |less

Thanks to Jed Brown for help developing this example.

