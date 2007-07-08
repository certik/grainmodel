all:
	./setpaths
	cd solver-petsc && make
	cd mesh && make
	cd assembly && make -j2

clean:
	cd solver-petsc && make distclean
	cd mesh && make clean
	cd assembly && make clean
	rm paths
	py.cleanup
