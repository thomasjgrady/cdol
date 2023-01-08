CC      = cc
CFLAGS  = -std=c99 -Wall -g3
SOURCES = libcdol.c
OBJECTS = $(SOURCES:.c=.o)
INCLUDE = -I .
LDFLAGS = -lc
TARGETS = libcdol_test.bin libcdol.so

all: $(TARGETS)

libcdol.so: $(SOURCES)	
	$(CC) $(CFLAGS) $(INCLUDE) -fPIC -shared -o $@ $(SOURCES) $(LDFLAGS)

libcdol_test.bin: $(SOURCES) libcdol.so
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ libcdol_test.c -L. -lcdol

.PHONY: clean install uninstall

clean:
	rm -v $(TARGETS)

PREFIX ?= /usr/local/lib

install:
	cp libcdol.so $(PREFIX)/libcdol.so

uninstall:
	rm -v $(PREFIX)/libcdol.so
