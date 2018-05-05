FC90 = mpif90 -fopenmp -fbounds-check
vector_files = mod_file.f90 main.f90

all:main_exe

main_exe:$(vector_files)
	$(FC90)	$(vector_files)	-o	$@
clean:
	rm	*_exe	*.mod	*~
