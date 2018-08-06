FILES = src/*.cpp lib/*.cpp lib/aux/*.cpp
all:
	mkdir -p bin
	mpic++ -DMPI -std=c++14 ${FILES} -o bin/loki
clean:
	rm bin/loki