CFLAGS=-c -O3 -I$(BOOST_INCLUDEDIR) -std=gnu++11

all: generate infer

generate: generate.o
	g++ -L$(LD_LIBRARY_PATH) generate.o -o generate

generate.o: generate.cpp Particle.hpp
	g++ $(CFLAGS) generate.cpp

infer: infer.o SmoothKernelApproximation.o
	g++ -L$(LD_LIBRARY_PATH) infer.o SmoothKernelApproximation.o -o infer

infer.o: infer.cpp Particle.hpp
	g++ $(CFLAGS) infer.cpp
	
SmoothKernelApproximation.o: SmoothKernelApproximation.cpp Particle.hpp
	g++ $(CFLAGS) SmoothKernelApproximation.cpp

load:
	module load boost
	module load gcc/4.8.3

clean:
	rm *.o generate infer
