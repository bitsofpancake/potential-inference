#include <iostream>
#include <cmath>
#include <random>
#include <boost/numeric/odeint.hpp>
#include "SmoothKernelApproximation.hpp"
#include "Particle.hpp"

using namespace boost::numeric::odeint;
/*
double loglikelihood(const std::vector<Particle> &data, const HamiltonianSystem &sys) {
	// Generate the distribution function
	const double T = 20.0;
	SmoothKernelApproximation f;
	f.add(data);
	
	const int samplesPerUnitTime = 10;
	const double dt = 1.0 / samplesPerUnitTime;
	symplectic_rkn_sb3a_mclachlan<vector_t> stepper;
	std::vector<Particle> evolved_data = data;
	int i = 0;
	int n = data.size();
	for (double t = 0.0; t < T; t += dt) {
		#pragma omp parallel for
		for (int i = 0; i < n; i++)
			stepper.do_step(sys, evolved_data[i].q, evolved_data[i].p, t, dt);
	
		if (++i % samplesPerUnitTime == 0)
			f.add(evolved_data);
	}
	f.save();
	
	// Multiply likelihoods together.
	double loglikelihood = 0.0;
	#pragma omp parallel for reduction(+:loglikelihood)
	for (int i = 0; i < n; i++)
		loglikelihood += log(f(data[i]));
	return loglikelihood;
}
*/
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
		
		// Advance t -- necessary because of the multithreading. 
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
		loglikelihood += log(likelihoods[i]);
	return loglikelihood;
}

double annealing_schedule(int t) {
	return 1.0 * exp(-t / 500.0);
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
	int iteration = 0;
	while (true) {
		param_t candidate_param;
		std::normal_distribution<double> jumping(0, annealing_schedule(iteration));
		for (int i = 0; i < candidate_param.size(); i++)
			candidate_param[i] = jumping(engine) + current_param[i];
	
		double candidate_loglikelihood = log(prior(candidate_param)) + loglikelihood(data, HamiltonianSystem(candidate_param));
		double acceptance_ratio = exp(candidate_loglikelihood - current_loglikelihood);
		if (acceptance_ratio >= acceptance(engine)) {
			//std::cout << candidate_param << std::endl;
			current_param = candidate_param;
			current_loglikelihood = candidate_loglikelihood;
			iteration++;
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