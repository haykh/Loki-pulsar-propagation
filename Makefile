FILES = src/*.cpp lib/*.cpp lib/aux/*.cpp

CC = h5c++
CFLAGS = -std=c++14

FLAGS = -DINTBACK

all:
	mkdir -p bin
	${CC} ${CFLAGS} ${FLAGS} ${FILES} -o bin/loki

clean:
	rm bin/loki
