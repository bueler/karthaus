ALL: ssaflowline

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables

ssaflowline : ssaflowline.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB}
	rm -f ssaflowline.o

# version zero is for pedagogical purposes; it shows the first-running-state
#   of the PETSc implementation
ssaflowline_v0 : ssaflowline_v0.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB}
	rm -f ssaflowline_v0.o

# these tests use a fine grid but run fast anyway
test: ssaflowline
	./ssaflowline -da_grid_x 20000 -ssa_fd -mat_fd_type ds -snes_max_it 500
	./ssaflowline -da_grid_x 20000 -ssa_fd -mat_fd_type ds -snes_max_it 500 \
	   -ssa_guess 2
	mpiexec -n 4 ./ssaflowline -da_grid_x 20000
	mpiexec -n 4 ./ssaflowline -da_grid_x 20000 -snes_mf_operator \
	    -mat_mffd_type ds -snes_ls basic -ssa_picard

