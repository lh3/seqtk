CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function

all:seqtk tabtk

seqtk:seqtk.c khash.h kseq.h
		$(CC) $(CFLAGS) seqtk.c -o $@ -lz -lm

tabtk:tabtk.c kseq.h kstring.h
		$(CC) $(CFLAGS) tabtk.c -o $@ -lz -lm

clean:
		rm -fr gmon.out *.o ext/*.o a.out seqtk tabtk *~ *.a *.dSYM session*
