#include <iostream>
#include <cmath>
#include <random>
#include <boost/numeric/odeint.hpp>
#include "SmoothKernelApproximation.hpp"
#include "Particle.hpp"
#include "AdaptiveMetropolisHastings.hpp"

using namespace boost::numeric::odeint;

double loglikelihood(const std::vector<Particle> &data, const HamiltonianSystem &sys) {
	// Generate the distribution function
	const double T = 20.0;
	const double dt = 0.1;
	const double sampleRate = 1;
	int n = data.size();
	
	std::vector<Particle> evolved_data = data;
	std::vector<double> likelihoods(n);
	SmoothKernelApproximation f;
	f.add(evolved_data);
	f.save();
	
	#pragma omp parallel for
	for (int i = 0; i < n; i++)
		likelihoods[i] = f(data[i]);
	
	symplectic_rkn_sb3a_mclachlan<vector_t> stepper;
	double t = 0;
	while (t <= T) {
		double tf = t + 1.0 / sampleRate;
		
		#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			double tt = t;
			while (tt <= tf) {
				stepper.do_step(sys, evolved_data[i].q, evolved_data[i].p, tt, dt);
				tt += dt;
			}
		}
		
		// Advance t. 
		while (t <= tf)
			t += dt;
		
		// Add the evolved data.
		SmoothKernelApproximation f;
		f.add(evolved_data);
		f.save();
		
		#pragma omp parallel for
		for (int i = 0; i < n; i++)
			likelihoods[i] += f(data[i]);
	}
	
	// Multiply likelihoods together.
	double loglikelihood = 0.0;
	for (int i = 0; i < n; i++)
		loglikelihood += log(likelihoods[i] / T);
	return loglikelihood;
}

/*
// Use this to benchmark
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
	
	param_t current_param;
	for (int i = 0; i < current_param.size(); i++)
		current_param[i] = argc - 1 > i ? atof(argv[1 + i]) : 1;
	std::cout << loglikelihood(data, HamiltonianSystem(current_param)) << std::endl;
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

//*
	// Metropolis-Hastings
	param_t current_param;
	
	// Initialize arbitrary starting point.
	for (int i = 0; i < current_param.size(); i++)
		current_param[i] = argc - 1 > i ? atof(argv[1 + i]) : 1;
			
	std::random_device rd;
	std::mt19937 engine(rd());
	std::uniform_real_distribution<double> acceptance(0, 1);
	double current_loglikelihood = log(prior(current_param)) + loglikelihood(data, HamiltonianSystem(current_param));
	
	AdaptiveMetropolisHastingsProposalDistribution<param_size> proposal;
	proposal.add(current_param.data());
	while (true) {
		param_t candidate_param;
		proposal.next(current_param.data(), candidate_param.data());
	
		double candidate_loglikelihood = log(prior(candidate_param)) + loglikelihood(data, HamiltonianSystem(candidate_param));
		double acceptance_ratio = exp(candidate_loglikelihood - current_loglikelihood);
		//std::cerr << acceptance_ratio << "\te^" << (candidate_loglikelihood - current_loglikelihood) << std::endl;
		if (acceptance_ratio >= acceptance(engine)) {
			//std::cout << candidate_param << std::endl;
			current_param = candidate_param;
			current_loglikelihood = candidate_loglikelihood;
			std::cerr << "Current loglikelihood: " << current_loglikelihood << "; param: " << current_param << std::endl;
			proposal.add(current_param.data());
		} else {
			//std::cout << candidate_param << " (rejected)" << std::endl;
		}
		std::cout << candidate_param << candidate_loglikelihood << std::endl;
	}
	
/*/	
	// Calculate likelihoods for all alpha.
	for (double alpha = 1.0; alpha < 3.0; alpha += 0.1) {
		param_t current_param;
		current_param[0] = alpha;
		std::cout << alpha << '\t' << loglikelihood(data, HamiltonianSystem(current_param)) << std::endl;
	}
//*/
	
	//std::cout << loglikelihood(data, HamiltonianSystem(1.8)) << std::endl;
}