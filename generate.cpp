#include <iostream>
#include <cmath>
#include <random>
#include <boost/numeric/odeint.hpp>
#include "Particle.hpp"

using namespace boost::numeric::odeint;

void normalize_vector(vector_t &v) {
	double r = 0;
	for (int i = 0; i < dim; i++)
		r += v[i]*v[i];
	r = sqrt(r);
	for (int i = 0; i < dim; i++)
		v[i] /= r;
}

int main(int argc, char *argv[]) {
	const int n = atoi(argv[1]);
	param_t param;
	for (int i = 0; i < param.size(); i++)
		param[i] = atof(argv[2 + i]);
	HamiltonianSystem sys(param);
	
	const double dt = 0.1;
	const double mixingTime = 10000.0; // When in doubt, increase this!
	
	// Generate initial conditions by randomly populating stars, and then evolving them until phase-mixed.
	#pragma omp parallel
	{
		std::random_device rd;
		std::mt19937 engine(rd());
		std::uniform_real_distribution<double> rand(0.25, 1.5);
		std::discrete_distribution<int> rand_sign {-1, 1};
		
		#pragma omp for
		for (int j = 0; j < n; j++) {
			Particle particle;
			
			for (int i = 0; i < 2*dim; i++)
				particle.coords[i] = rand_sign(engine) * rand(engine);
			
			symplectic_rkn_sb3a_mclachlan<vector_t> stepper;
			for (double t = 0.0; t <= mixingTime; t += dt)
				stepper.do_step(sys, particle.q, particle.p, t, dt);
			
			#pragma omp critical
			{
				for (int i = 0; i < 2*dim; i++)
					std::cout << particle.coords[i] << '\t';
				std::cout << std::endl;
			}
		}
	}
}