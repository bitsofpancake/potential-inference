#include <iostream>
#include <cmath>
#include <random>
#include <boost/numeric/odeint.hpp>
#include "SmoothKernelApproximation.hpp"
#include "Particle.hpp"

using namespace boost::numeric::odeint;

double loglikelihood(std::vector<Particle> data, HamiltonianSystem sys) {
	// Generate the distribution function
	const double T = 1000.0;
	SmoothKernelApproximation f;
	f.add(data);
	
	const int samplesPerUnitTime = 10;
	const double dt = 1.0 / samplesPerUnitTime;
	symplectic_rkn_sb3a_mclachlan<vector_t> stepper;
	int i = 0;
	int n = data.size();
	for (double t = 0.0; t < T; t += dt) {
		#pragma omp parallel for
		for (int i = 0; i < n; i++)
			stepper.do_step(sys, data[i].q, data[i].p, t, dt);
	
		if (++i % samplesPerUnitTime == 0)
			f.add(data);
	}
	f.save();
	
	// Multiply likelihoods together.
	double loglikelihood = 0.0;
	#pragma omp parallel for reduction(+:loglikelihood)
	for (int i = 0; i < n; i++)
		loglikelihood += log(f(data[i]));
	return loglikelihood / T;
}

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

	// Metropolis-Hastings
	param_t current_param;
	
	// Initialize arbitrary starting point.
	for (int i = 0; i < current_param.size(); i++)
		current_param[i] = 1;
	
	std::random_device rd;
	std::mt19937 engine(rd());
	std::uniform_real_distribution<double> acceptance(0, 1);
	std::normal_distribution<double> jumping(0, 0.1);
	double current_loglikelihood = log(prior(current_param)) + loglikelihood(data, HamiltonianSystem(current_param));
	while (true) {
		param_t candidate_param;
		for (int i = 0; i < candidate_param.size(); i++)
			candidate_param[i] = jumping(engine) + current_param[i];
	
		double candidate_loglikelihood = log(prior(candidate_param)) + loglikelihood(data, HamiltonianSystem(candidate_param));
		double acceptance_ratio = exp(candidate_loglikelihood - current_loglikelihood);
		if (acceptance_ratio >= acceptance(engine)) {
			//std::cout << candidate_param << std::endl;
			current_param = candidate_param;
			current_loglikelihood = candidate_loglikelihood;
		} else {
			//std::cout << candidate_param << " (rejected)" << std::endl;
		}
		std::cout << candidate_param << candidate_loglikelihood << std::endl;
	}
/*/	
	// Calculate likelihoods for all alpha.
	param_t current_param;
	for (double alpha = 1.0; alpha < 3.0; alpha += 0.1) {
		current_param[0] = alpha;
		std::cout << alpha << '\t' << loglikelihood(data, HamiltonianSystem(current_param)) << std::endl;
	}
*/
	
	//std::cout << loglikelihood(data, HamiltonianSystem(1.8)) << std::endl;
}