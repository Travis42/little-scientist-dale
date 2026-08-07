# DALE: Dual-seed Algorithm for Latent Enumeration
# Pure C99, only libm. `make` builds ./benchmark.

CC      ?= cc
CFLAGS  ?= -std=c99 -O2 -Wall -Wextra -pedantic
LDFLAGS ?=
LDLIBS  ?= -lm

OBJS = dale.o benchmark.o

all: benchmark

benchmark: $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

dale.o:    dale.c dale.h
	$(CC) $(CFLAGS) -c dale.c -o $@

benchmark.o: benchmark.c dale.h
	$(CC) $(CFLAGS) -c benchmark.c -o $@

check: benchmark
	./benchmark --data ../example --tf SP1

clean:
	rm -f benchmark $(OBJS)

.PHONY: all check clean
