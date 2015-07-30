#include <iostream>
#include <cmath>
#include <random>
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

double sign(const double x) {
	return std::signbit(x) ? -1.0 : 1.0;
}

const int dim = 1;
typedef boost::array<double, dim> vector_t;
typedef vector_t particle_t[2];

struct Particle {
	vector_t p;
	vector_t q;
};

class HamiltonianSystem {
	const double alpha;

	public:
		HamiltonianSystem(const double alpha) : alpha(alpha) {};

		// Assume the position derivative is p.
		// Returns the momentum derivative.
		void operator()(const vector_t &q, vector_t &dpdt) const {
			dpdt[0] = -0.5 * pow(fabs(q[0]), alpha - 1) * sign(q[0]);
		}
};

class SmoothKernelApproximation {
	std::vector<Particle> data;
	
	double kernel() {
	
	}
	
	public:
		void add(const std::vector<Particle> others) {
			data.insert(data.end(), others.begin(), others.end());
		};
		
		double operator()(Particle particle) const {
			return this->operator()(particle.p, particle.q);
		};
		
		double operator()(vector_t p, vector_t q) const {
			return 1.0 / data.size();
		};
};

double loglikelihood(std::vector<Particle> data, HamiltonianSystem sys) {
	// Generate the distribution function
	const int T = 20;
	SmoothKernelApproximation f;
	f.add(data);
	
	const double dt = 0.1;
	symplectic_rkn_sb3a_mclachlan<vector_t> stepper;
	for (double t = 0.0; t < T; t += dt) {
		for (Particle &particle : data)
			stepper.do_step(sys, particle.q, particle.p, t, dt);
		
		f.add(data);
	}
	
	// Multiply likelihoods together.
	double loglikelihood = 0.0;
	for (Particle &particle : data)
		loglikelihood += log(f(particle));
	return loglikelihood;
}

int main(int argc, char *argv[]) {
	HamiltonianSystem sys(1.5);
	const int n = atoi(argv[1]);
	const double dt = 0.1;
	
	// Generate initial conditions.
	std::uniform_real_distribution<double> rand(-15, 15);
	std::default_random_engine engine;
	std::vector<Particle> data(n);
	for (int i = 0; i < n; i++) {
		data[i].q[0] = rand(engine);
		data[i].p[0] = rand(engine);
		
		symplectic_rkn_sb3a_mclachlan<vector_t> stepper;
		for (double t = 0.0; t < 100; t += dt)
			stepper.do_step(sys, data[i].q, data[i].p, t, dt);
	}
	
	std::cout << loglikelihood(data, HamiltonianSystem(1.6));
}