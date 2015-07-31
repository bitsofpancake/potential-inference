CFLAGS=-c -O3 -I$(BOOST_INCLUDEDIR) -std=gnu++11

all: main.o SmoothKernelApproximation.o
	g++ -L$(LD_LIBRARY_PATH) main.o SmoothKernelApproximation.o -o main

main.o: main.cpp
	g++ $(CFLAGS) main.cpp
	
SmoothKernelApproximation.o: SmoothKernelApproximation.cpp
	g++ $(CFLAGS) SmoothKernelApproximation.cpp

load:
	module load boost
	module load gcc/4.8.3

clean:
	rm *.o main
