#ifndef PARTICLE_HPP
#define PARTICLE_HPP

const int dim = 1;
typedef double vector_t[dim];
union Particle {
	double coords[dim*2];
	struct {
		vector_t p;
		vector_t q;
	};
};

#endif