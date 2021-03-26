CC ?= gcc
CFLAGS ?= -g -Wall -O2 -Wno-unused-function
PREFIX ?= /usr/local
BINDIR ?= $(DESTDIR)$(PREFIX)/bin
INSTALL ?= install
MKDIR ?= mkdir

all:seqtk

seqtk:seqtk.c khash.h kseq.h
		$(CC) $(CFLAGS) seqtk.c -o $@ -lz -lm

install:all
		$(MKDIR) -p $(BINDIR)
		$(INSTALL) seqtk $(BINDIR)

clean:
		rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session*
