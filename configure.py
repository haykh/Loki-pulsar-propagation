import sys

arguments = sys.argv[1:]
makefile = "FILES = src/*.cpp lib/*.cpp lib/aux/*.cpp\n"
makefile += "all:\n\tmkdir -p bin\n" 

flags = ""
cppflags = "-std=c++14 "

if "-intback" in arguments:
	flags += "-DINTBACK "
if "-findnan" in arguments:
	flags += "-DFINDNAN "

if "-mpi" in arguments:
	makefile += "\tmpic++ -DMPI "
else:
	makefile += "\tg++ "

makefile += cppflags + flags + "${FILES} -o bin/loki\n"
#makefile += "\trm *.o\n\n"

makefile += "clean:\n\trm bin/loki"

text_file = open("Makefile", "w")
text_file.write(makefile)
text_file.close()
