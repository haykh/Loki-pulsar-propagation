FILES = src/*.cpp lib/*.cpp lib/aux/*.cpp

CC = g++
CFLAGS = -std=c++14
MPICC = mpic++ -DMPI

FLAGS =# -DINTBACK

all:
	mkdir -p bin
	${CC} ${CFLAGS} ${FLAGS} ${FILES} -o bin/loki

clean:
	rm bin/loki

mpi:
	mkdir -p bin
	${MPICC} ${CFLAGS} ${FLAGS} ${FILES} -o bin/loki
