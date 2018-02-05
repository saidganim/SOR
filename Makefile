
all: sor-seq sor-par

run: all
	par/sor-par ${NSIZE} > par/output
	seq/sor-seq ${NSIZE} > seq/output

diff: run
	diff -s par/output seq/output > diff_par_seq
	diff -y par/output seq/output >> diff_par_seq
sor-seq: seq/
	$(MAKE) -C seq

sor-par: par/
	$(MAKE) -C par

clean:
	rm par/sor-par seq/sor-seq
	rm par/output* seq/output*
	rm diff_*
