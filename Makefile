

MPICC=mpic++
MPIRUN=mpirun

run: default
	$(MPIRUN) bin.elf

default: sor-par.cpp
	$(MPICC) sor-par.cpp -o bin.elf

clean: bin.elf
	rm bin.elf
