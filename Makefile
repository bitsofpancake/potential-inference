CFLAGS=-c -O3 -funroll-loops -fopenmp -I$(BOOST_INCLUDEDIR) -std=gnu++11
LINKERFLAGS=-fopenmp

all: generate infer

generate: generate.o
	g++ -L$(LD_LIBRARY_PATH) $(LINKERFLAGS) generate.o -o generate

generate.o: generate.cpp Particle.hpp
	g++ $(CFLAGS) generate.cpp

infer: infer.o SmoothKernelApproximation.o
	g++ -L$(LD_LIBRARY_PATH) $(LINKERFLAGS) infer.o SmoothKernelApproximation.o -o infer

infer.o: infer.cpp Particle.hpp AdaptiveMetropolisHastings.hpp
	g++ $(CFLAGS) infer.cpp
	
SmoothKernelApproximation.o: SmoothKernelApproximation.cpp Particle.hpp
	g++ $(CFLAGS) SmoothKernelApproximation.cpp

load:
	module load boost
	module load gcc/4.8.3

clean:
	rm -f *.o generate infer

bench: clean infer
	time ./infer < data200k

benchg: clean generate
	time ./generate 1000 1.8 > /dev/null