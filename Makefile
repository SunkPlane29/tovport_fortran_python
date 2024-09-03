.PHONY: runplot
runplot: runall plot

.PHONY: runall
runall:
	cd fortran && make run
	cd python && make run
	cd julia && make run

.PHONY: plot
plot: initout
	python3 plots.py

.PHONY: initout
initout:
	mkdir -p out
