# spec_f
Speciation scripts using fortran (and a little of Matlab)

## Description of relevant files
- 'utils.f95' contains user-defined functions
- 'main.f95' computes the model (it needs to compile 'utils.f95' first)
	```
	gfortran -c utils.f95
	gfortran -o main.out main.f95 utils.o
	./main.out
	```

- 'main_full.f95' computes the model (contains local functions)
	```
	gfortran -o main_full.out main_full.f95
	./main_full.out
	```

## How to plot?
Running 'main.f95' or 'main_full.f95' creates files 'pop.dat' (last population matrix) and 'vmat.dat' (mating segment length), both readable by matlab. Run 'grph_fortran.m' in MATLAB to get the graphs.



