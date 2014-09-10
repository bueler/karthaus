ALL: ssaflowline

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables

ssaflowline : ssaflowline.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB}
	rm -f ssaflowline.o

# these tests use a fine grid but run fast anyway
test: ssaflowline
	mpiexec -n 1 ./ssaflowline -da_grid_x 20001
	mpiexec -n 4 ./ssaflowline -da_grid_x 20001
	mpiexec -n 4 ./ssaflowline -da_grid_x 20001 -ssa_guess 2
	mpiexec -n 4 ./ssaflowline -da_grid_x 20001 -snes_mf_operator -ssa_picard

