

MPICC=mpic++
MPIRUN=mpirun

run: default
	$(MPIRUN) bin.elf

default: main.cpp
	$(MPICC) main.cpp -o bin.elf

clean: bin.elf
	rm bin.elf
