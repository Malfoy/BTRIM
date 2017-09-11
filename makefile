CC=g++
CFLAGS=  -Wall -Wextra  -Ofast -std=c++11  -pthread -pipe -fopenmp
LDFLAGS=-pthread -fopenmp


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O3  -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O3 -g
LDFLAGS=-g
endif


EXEC=btrim bubble badvisor

all: $(EXEC)

bubble:   bubble.o
	$(CC) -o $@ $^ $(LDFLAGS)

bubble.o:   btrimBubble.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

badvisor:   badvisor.o
	$(CC) -o $@ $^ $(LDFLAGS)

badvisor.o:   Badvisor.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

btrim:   btrim.o
	$(CC) -o $@ $^ $(LDFLAGS)

btrim.o:   btrim.cpp
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
