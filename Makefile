FILES = src/*.cpp lib/*.cpp lib/aux/*.cpp

CC = g++
CFLAGS = -std=c++14

FLAGS = -DINTBACK

all:
	mkdir -p bin
	${CC} ${CFLAGS} ${FLAGS} ${FILES} -o bin/loki

clean:
	rm bin/loki
