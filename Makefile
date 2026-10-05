CC = gcc
CFLAGS = -std=gnu99 -Wall -O2
LDLIBS = -lm

HPART_OBJS = hpart.o readpart.o greedy.o

.PHONY: all clean

all: hpart hmake

hpart: $(HPART_OBJS)
	$(CC) $(CFLAGS) -o $@ $(HPART_OBJS) $(LDLIBS)

hmake: hmake.o
	$(CC) $(CFLAGS) -o $@ hmake.o

hpart.o: hpart.c hpart.h
readpart.o: readpart.c hpart.h
greedy.o: greedy.c hpart.h
hmake.o: hmake.c hmake.h

clean:
	rm -f hpart hmake *.o
