FILES = src/*.cpp lib/*.cpp

CC = g++
CFLAGS = -std=c++14

MPICC = mpic++ -DMPI

all:
	mkdir -p bin
	${CC} ${CFLAGS} ${FILES} -o bin/loki

clean:
	rm bin/loki

mpi:
	mkdir -p bin
	${MPICC} ${CFLAGS} ${FILES} -o bin/loki
