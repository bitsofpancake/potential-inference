#include <iostream>
#include <cmath>
#include <random>
#include <boost/numeric/odeint.hpp>
#include "Particle.hpp"

using namespace boost::numeric::odeint;

double sign(const double x) {
	return std::signbit(x) ? -1.0 : 1.0;
}

typedef double param_t;
class HamiltonianSystem {
	const param_t alpha;

 public:
	HamiltonianSystem(const param_t alpha) : alpha(alpha) {};

	// Assume the position derivative is p.
	// Returns the momentum derivative.
	void operator()(const vector_t &q, vector_t &dpdt) const {
		dpdt[0] = -0.5 * pow(fabs(q[0]), alpha - 1) * sign(q[0]);
	};
};

int main(int argc, char *argv[]) {
	const int n = atoi(argv[1]);
	HamiltonianSystem sys(atof(argv[2]));
	const double dt = 0.1;
	const double mixingTime = 500.0; // When in doubt, increase this!
	
	// Generate initial conditions by randomly populating stars, and then evolving them until phase-mixed.
	std::random_device rd;
	std::mt19937 engine(rd());
	std::uniform_real_distribution<double> rand(5, 15);
	std::discrete_distribution<int> rand_sign {-1, 1};
	std::vector<Particle> data(n);
	for (Particle &p : data) {
		p.q[0] = rand_sign(engine) * rand(engine);
		p.p[0] = rand_sign(engine) * rand(engine);
		
		symplectic_rkn_sb3a_mclachlan<vector_t> stepper;
		for (double t = 0.0; t <= mixingTime; t += dt)
			stepper.do_step(sys, p.q, p.p, t, dt);
		
		for (int i = 0; i < 2*dim; i++)
			std::cout << p.coords[i] << '\t';
		std::cout << std::endl;
	}
}