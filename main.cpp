#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <boost/numeric/odeint.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

using namespace boost::accumulators;
using namespace boost::numeric::odeint;

double sign(const double x) {
	return std::signbit(x) ? -1.0 : 1.0;
}

const int dim = 1;
typedef double vector_t[dim];
struct Particle {
	vector_t p;
	vector_t q;
};
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

class SmoothKernelApproximation {
	std::vector<Particle> data;
	mutable accumulator_set<double, stats<tag::variance>> accumulator_p[dim];
	mutable accumulator_set<double, stats<tag::variance>> accumulator_q[dim];
	double bandwidth_p[dim];
	double bandwidth_q[dim];
	double vol = 1.0;
	
	double kernel(const double u) const {
		if (u <= -1 || u >= 1)
			return 0;
		return 0.75 * (1 - u*u);
	}
	
	public:
		void add(const std::vector<Particle> others) {
			data.insert(data.end(), others.begin(), others.end());
			for (const Particle &particle : others) {
				for (int i = 0; i < dim; i++) {
					accumulator_p[i](particle.p[i]);
					accumulator_q[i](particle.q[i]);
				}
			}
		};
		
		void save() {
			// Calculate bandwidth.
			for (int i = 0; i < dim; i++) {
				// Scott's rule
				bandwidth_p[i] = sqrt(variance(accumulator_p[i])) * pow(data.size(), -1.0 / (4 + 2 * dim));
				bandwidth_q[i] = sqrt(variance(accumulator_q[i])) * pow(data.size(), -1.0 / (4 + 2 * dim));
				vol *= bandwidth_p[i] * bandwidth_q[i];
			}
		}
		
		double operator()(const Particle a) const {
			double value = 0.0;			
			for (const Particle &b : data) {
				double u_sq = 0;
				for (int i = 0; i < dim; i++) {
					u_sq += pow((a.p[i] - b.p[i]) / bandwidth_p[i], 2);
					u_sq += pow((a.q[i] - b.q[i]) / bandwidth_q[i], 2);
				}
				value += kernel(sqrt(u_sq));
			}
				
			return value / data.size() / vol;
		};
};

double loglikelihood(std::vector<Particle> data, HamiltonianSystem sys) {
	// Generate the distribution function
	const double T = 20.0;
	SmoothKernelApproximation f;
	f.add(data);
	
	const double dt = 0.01;
	symplectic_rkn_sb3a_mclachlan<vector_t> stepper;
	int i = 0;
	for (double t = 0.0; t < T; t += dt) {
		for (Particle &particle : data)
			stepper.do_step(sys, particle.q, particle.p, t, dt);
		
		if (++i % 100 == 0)
			f.add(data);
	}
	f.save();
	
	// Multiply likelihoods together.
	double loglikelihood = 0.0;
	for (const Particle &particle : data)
		loglikelihood += log(f(particle));
	return loglikelihood;
}

int main(int argc, char *argv[]) {
	const int n = atoi(argv[1]);
	HamiltonianSystem sys(atof(argv[2]));
	const double dt = 0.01;
	
	// Generate initial conditions by randomly populating stars, and then evolving them until phase-mixed.
	std::random_device engine;
	std::uniform_real_distribution<double> rand(5, 15);
	std::discrete_distribution<int> rand_sign {-1, 1};
	std::vector<Particle> data(n);
	for (int i = 0; i < n; i++) {
		data[i].q[0] = rand_sign(engine) * rand(engine);
		data[i].p[0] = rand_sign(engine) * rand(engine);
		
		symplectic_rkn_sb3a_mclachlan<vector_t> stepper;
		for (double t = 0.0; t < 100; t += dt)
			stepper.do_step(sys, data[i].q, data[i].p, t, dt);
	}
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