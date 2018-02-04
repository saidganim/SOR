

MPICC=mpic++
MPIRUN=mpirun

run: sor-par
	$(MPIRUN) -np 4 sor-par 100

sor-par: sor-par.cpp
	$(MPICC) sor-par.cpp -o sor-par

clean: sor-par
	rm sor-par
