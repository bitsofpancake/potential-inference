#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <boost/numeric/odeint.hpp>
#include "SmoothKernelApproximation.hpp"
#include "Particle.hpp"

using namespace boost::numeric::odeint;

double sign(const double x) {
	return std::signbit(x) ? -1.0 : 1.0;
}

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

double loglikelihood(std::vector<Particle> data, HamiltonianSystem sys) {
	// Generate the distribution function
	const double T = 20.0;
	SmoothKernelApproximation2 f;
	f.add(data);
	
	const int samplesPerUnitTime = 10;
	const double dt = 1.0 / samplesPerUnitTime;
	symplectic_rkn_sb3a_mclachlan<vector_t> stepper;
	int i = 0;
	for (double t = 0.0; t < T; t += dt) {
		for (Particle &particle : data)
			stepper.do_step(sys, particle.q, particle.p, t, dt);
		
		if (++i % samplesPerUnitTime == 0)
			f.add(data);
	}
	f.save();
	
	// Multiply likelihoods together.
	double loglikelihood = 0.0;
	for (const Particle &particle : data)
		loglikelihood += log10(f(particle));
	return loglikelihood;
}
/*
int main(int argc, char *argv[]) {
	const int n = atoi(argv[1]);
	HamiltonianSystem sys(atof(argv[2]));
	const double dt = 0.1;
	const double mixingTime = 500.0; // When in doubt, increase this!
	
	// Generate initial conditions by randomly populating stars, and then evolving them until phase-mixed.
	std::random_device engine;
	std::uniform_real_distribution<double> rand(5, 15);
	std::discrete_distribution<int> rand_sign {-1, 1};
	std::vector<Particle> data(n);
	for (int i = 0; i < n; i++) {
		data[i].q[0] = rand_sign(engine) * rand(engine);
		data[i].p[0] = rand_sign(engine) * rand(engine);
		
		symplectic_rkn_sb3a_mclachlan<vector_t> stepper;
		for (double t = 0.0; t <= mixingTime; t += dt)
			stepper.do_step(sys, data[i].q, data[i].p, t, dt);
		
		std::cout << data[i].q[0] << '\t' << data[i].p[0] << std::endl;
	}
}
*/

int main(int argc, char *argv[]) {

	// Read data from stdin
	std::vector<Particle> data;
	while (true) {
		Particle particle;
		bool success = false;
		for (int i = 0; i < 2*dim; i++)
			success = std::cin >> particle.coords[i];
		if (!success)
			break;
		data.push_back(particle);
	}
	/*
	SmoothKernelApproximation f;
	f.add(data);
	f.save();
	f(data[0]);
	
	std::cout << std::endl;
	
	SmoothKernelApproximation2 f2;
	f2.add(data);
	f2.save();
	f2(data[0]);
	*/
	/*
	
	for (int i = 0; i < data.size(); i++) {
		std::cout << i << '\t' << f(data[i]) << '\t' << f2(data[i]) << std::endl;
		
	}*/
	
	/*
	// Metropolis-Hastings
	param_t current_param = 2.0;
	std::uniform_real_distribution<double> acceptance(0, 1);
	double current_loglikelihood = loglikelihood(data, HamiltonianSystem(current_param));
	while (true) {
		std::normal_distribution<param_t> jumping(current_param, 0.1);
		param_t candidate_param = jumping(engine);
	
		double candidate_loglikelihood = loglikelihood(data, HamiltonianSystem(candidate_param));
		double acceptance_ratio = exp(candidate_loglikelihood - current_loglikelihood);
		if (acceptance_ratio >= acceptance(engine)) {
			std::cout << candidate_param << std::endl;
			current_param = candidate_param;
			current_loglikelihood = candidate_loglikelihood;
		} else {
			//std::cout << candidate_param << " (rejected)" << std::endl;
		}
	}
	*/
	
	// Calculate likelihoods for all alpha.
	for (double alpha = 1.0; alpha < 3.0; alpha += 0.1)
		std::cout << alpha << '\t' << loglikelihood(data, HamiltonianSystem(alpha)) << std::endl;
}