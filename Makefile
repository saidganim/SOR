
all: sor-seq sor-par

run: all
	par/sor-par ${NSIZE} > par/output
	seq/sor-seq ${NSIZE} > seq/output

sor-seq: seq/
	$(MAKE) -C seq

sor-par: par/
	$(MAKE) -C par

clean:
	rm par/sor-par
	rm seq/sor-seq
