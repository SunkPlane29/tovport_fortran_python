.PHONY: initout
initout:
	mkdir -p out

.PHONY: build
build: initout
	gfortran -I /usr/local/lib/lapack95/lapack95_modules -c tov.f90 spline.f90 diff.f90 constants.f90
	gfortran -L /usr/local/lib/lapack95 -L /usr/lib/x86_64-linux-gnu -L /usr/local/lib/lapack-3.10.1 \
		-o main.out main.f90 tov.o spline.o diff.o constants.o -l:lapack95.a -llapack -lblas -l:libtmglib.a
	rm -f *.o

.PHONY: run
run: build 
	./main.out
	