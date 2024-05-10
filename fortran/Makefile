build:
	gfortran -I /usr/local/lib/lapack95/lapack95_modules -c spline.f90
	gfortran -L /usr/local/lib/lapack95 -L /usr/lib/x86_64-linux-gnu -L /usr/local/lib/lapack-3.10.1 \
		-o main main.f90 spline.o -l:lapack95.a -llapack -lblas -l:libtmglib.a

run: build 
	./main
	