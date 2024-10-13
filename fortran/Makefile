NTHREADS := $(shell nproc --all)

.PHONY: initout
initout:
	mkdir -p out

.PHONY: build
build: initout
	gfortran -fopenmp -O3 -ftree-parallelize-loops=$(NTHREADS) -I ../dist/LAPACK95/lapack95_modules -c spline.f90 diff.f90 constants.f90
	gfortran -fopenmp -O3 -ftree-parallelize-loops=$(NTHREADS) -I ../dist/LAPACK95/lapack95_modules -c tov.f90
	gfortran -fopenmp -O3 -ftree-parallelize-loops=$(NTHREADS) -L ../dist/LAPACK95 -L /usr/lib/x86_64-linux-gnu -L /usr/lib/x86_64-linux-gnu/lapack -L /usr/lib/x86_64-linux-gnu/blas \
		-o main.out main.f90 tov.o spline.o diff.o constants.o -l:lapack95.a -llapack -lblas
	rm -f *.o *.mod

.PHONY: run
run: build
	time -p ./main.out
