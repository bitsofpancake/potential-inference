#ifndef PARTICLE_HPP
#define PARTICLE_HPP

const int dim = 2;
typedef double vector_t[dim];
union Particle {
	double coords[dim*2];
	struct {
		vector_t p;
		vector_t q;
	};
};

namespace { // Hacky, sorry

// Specify the potential here.
double sign(const double x) {
	return std::signbit(x) ? -1.0 : 1.0;
}

const int param_size = 2;
typedef std::array<double, param_size> param_t; // (R^2, q^2)
class HamiltonianSystem {
	const param_t param;

 public:
	HamiltonianSystem(const param_t param) : param(param) {};

	// Assume the position derivative is p.
	// Returns the momentum derivative.
	void operator()(const vector_t &q, vector_t &dpdt) const {
		
		// Logarithmic potential
		dpdt[0] = -q[0] * (param[0] + q[0]*q[0] + q[1]*q[1]/param[1]);
		dpdt[1] = -q[1]/param[1] * (param[0] + q[0]*q[0] + q[1]*q[1]/param[1]);
	
	//	dpdt[0] = -param[0] * pow(fabs(q[0]), param[0] - 1) * sign(q[0]);
	/*	
		// Toy galaxy potential
		double r_cubed_inv = pow(q[0]*q[0] + q[1]*q[1] + q[2]*q[2], -1.5);
		dpdt[0] = -q[0] * (param[0] + param[1] * r_cubed_inv);
		dpdt[1] = -q[1] * (param[0] + param[1] * r_cubed_inv);
		dpdt[2] = -q[2] * (param[0] + param[1] * r_cubed_inv);
	*/	
		// param[0] == mass of dark matter halo
		// param[1] == mass of black hole
	};
};

double prior(const param_t &param) {
	for (double p : param)
		if (p <= 0)
			return 0;
	/*if (param[1] <= 0.5 || param[1] >= 1.08*1.08)
		return 0;*/
	return 1;
};

std::ostream &operator<<(std::ostream &os, const param_t &param) {
	for (double p : param)
		os << p << '\t';
    return os;
};

}
#endif