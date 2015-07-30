CFLAGS=-c -I$(BOOST_INCLUDEDIR) -std=gnu++11

all: hello.o
	g++ -L$(LD_LIBRARY_PATH) hello.o -o hello

hello.o: hello.cpp
	g++ $(CFLAGS) hello.cpp
	
clean:
	rm *.o hello