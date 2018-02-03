

MPICC=mpic++
MPIRUN=mpirun

run: default
	$(MPIRUN) bin.elf

sor-par: sor-par.cpp
	$(MPICC) sor-par.cpp -o sor-par

clean: bin.elf
	rm bin.elf
