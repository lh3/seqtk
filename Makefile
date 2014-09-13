CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function

all:seqtk trimadap

seqtk:seqtk.c khash.h kseq.h
		$(CC) $(CFLAGS) seqtk.c -o $@ -lz -lm

trimadap:trimadap.c kseq.h ksw.h
		$(CC) $(CFLAGS) ksw.c trimadap.c -o $@ -lz -lm

clean:
		rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session*
