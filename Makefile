CC ?= gcc
CFLAGS+=-g -Wall -O2 -Wno-unused-function
BINDIR=/usr/local/bin

all: seqtk

seqtk: seqtk.c khash.h kseq.h
		$(CC) $(LDFLAGS) $(CFLAGS) $(CPPFLAGS) seqtk.c -o $@ -lz -lm

install: all
		install seqtk $(BINDIR)

clean:
		rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session*
