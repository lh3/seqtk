CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function

all:seqtk

version.h:
		echo "#define SEQTKVERSION \""$$(git describe --tags --always --dirty)"\"" > $@

seqtk:seqtk.c khash.h kseq.h version.h
		$(CC) $(CFLAGS) seqtk.c -o $@ -lz -lm
		rm version.h

clean:
		rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session*
