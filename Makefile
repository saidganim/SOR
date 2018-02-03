

MPICC=mpic++
MPIRUN=mpirun

run: sor-par
	$(MPIRUN) sor-par

sor-par: sor-par.cpp
	$(MPICC) sor-par.cpp -o sor-par

clean: sor-par
	rm sor-par
