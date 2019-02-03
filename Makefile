CC      = gcc
CFLAGS  = -g -Wall -O2 -Wno-unused-function
BINDIR  = /usr/local/bin

all:seqtk

seqtk:seqtk.c khash.h kseq.h constants.h seqtk.h
		$(CC) $(CFLAGS) seqtk.c -o $@ -lz -lm

install:all
		install seqtk $(BINDIR) $(INCL) $(SRC)

clean:
		rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session*
