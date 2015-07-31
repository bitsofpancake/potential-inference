CFLAGS=-c -O3 -I$(BOOST_INCLUDEDIR) -std=gnu++11

all: hello.o
	g++ -L$(LD_LIBRARY_PATH) hello.o -o hello

hello.o: hello.cpp
	g++ $(CFLAGS) hello.cpp

load:
	module load boost
	module load gcc/4.8.3

clean:
	rm *.o hello
