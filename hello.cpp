#include <iostream>
#include <cmath>
#include <random>
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

double dist(const Particle a, const Particle b) {
	double dist_sq = 0;
	for (int i = 0; i < 2 * dim; i++)
		dist_sq += pow(a.p[0] - b.p[0], 2) + pow(a.q[0] - b.q[0], 2);
	return sqrt(dist_sq);
}

class HamiltonianSystem {
	const double alpha;

	public:
		HamiltonianSystem(const double alpha) : alpha(alpha) {};

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
		
		double operator()(const Particle a) const {
			double value = 0.0;
			
			// Calculate bandwidth.
			double bandwidth_p[dim];
			double bandwidth_q[dim];
			double vol = 1.0;
			for (int i = 0; i < dim; i++) {
				bandwidth_p[i] = sqrt(variance(accumulator_p[i])) * pow(data.size(), -1.0 / (4 + 2 * dim));
				bandwidth_q[i] = sqrt(variance(accumulator_q[i])) * pow(data.size(), -1.0 / (4 + 2 * dim));
				vol *= bandwidth_p[i] * bandwidth_q[i];
			}
			
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
	for (const Particle &particle : data)
		loglikelihood += log(f(particle));
	return loglikelihood;
}

int main(int argc, char *argv[]) {
	const int n = atoi(argv[1]);
	HamiltonianSystem sys(atof(argv[2]));
	const double dt = 0.1;
	
	// Generate initial conditions by randomly populating stars, and then evolving them until phase-mixed.
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
	
	// Calculate likelihoods for all alpha.
	for (double alpha = 1.0; alpha < 3.0; alpha += 0.1)
		std::cout << alpha << '\t' << loglikelihood(data, HamiltonianSystem(alpha)) << std::endl;
}