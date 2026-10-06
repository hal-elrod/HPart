CC = gcc
CFLAGS = -std=gnu99 -Wall -O2
# Only hpart.c uses OpenMP (the parallel attempt loop in main()), so
# this is scoped to just its compile and hpart's final link step rather
# than added to CFLAGS for every file.
OMPFLAGS = -fopenmp
LDLIBS = -lm

HPART_OBJS = hpart.o readpart.o greedy.o

.PHONY: all clean

all: hpart hmake

hpart: $(HPART_OBJS)
	$(CC) $(CFLAGS) $(OMPFLAGS) -o $@ $(HPART_OBJS) $(LDLIBS) -lpthread

hmake: hmake.o
	$(CC) $(CFLAGS) -o $@ hmake.o

hpart.o: hpart.c hpart.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -c -o $@ hpart.c
readpart.o: readpart.c hpart.h
greedy.o: greedy.c hpart.h
hmake.o: hmake.c hmake.h

clean:
	rm -f hpart hmake *.o
