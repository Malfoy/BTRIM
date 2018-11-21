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


EXEC=btrim badvisor

all: $(EXEC)

btt:   btt.o
	$(CC) -o $@ $^ $(LDFLAGS)

btt.o:   btt.cpp
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
