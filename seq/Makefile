CFLAGS	:=
CFLAGS	+= -std=gnu99
CFLAGS	+= -Wall
CFLAGS	+= -O3
LDLIBS	:= -lm
BIN=./sor-seq

sor-seq: sor-seq.c
	gcc $(CFLAGS) -o sor-seq sor-seq.c $(LDLIBS)

clean:
	-rm -f *.o sor-seq *~ *core
