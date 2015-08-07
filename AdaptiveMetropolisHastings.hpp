#ifndef ADAPTIVE_METROPOLIS_HASTINGS_PROPOSAL_DISTRIBUTION_HPP
#define ADAPTIVE_METROPOLIS_HASTINGS_PROPOSAL_DISTRIBUTION_HPP

#include "Eigen/Dense"

// https://projecteuclid.org/download/pdf_1/euclid.bj/1080222083
template <int dim>
class AdaptiveMetropolisHastingsProposalDistribution {
	typedef Eigen::Matrix<double, dim, dim> Mat;
	typedef Eigen::Matrix<double, dim, 1> Vec;
	
	const double scale = 2.4 * 2.4 / dim; // Gelman, et al.
	const double epsilon = 0.01;
	const int initial_period = 3;
	
	Mat covariance;
	Vec sum;
	int samples = 0;
	
	mutable std::random_device rd;
	mutable std::mt19937 engine;
	mutable std::normal_distribution<double> normal_sample;
	Vec multivariate_normal_sample(const Mat &cov) const {
		Vec random;
		for (int i = 0; i < dim; i++)
			random[i] = normal_sample(engine);
		
		return cov.llt().matrixL() * random;
	}
	
	inline Mat outer(const Vec &v) const {
		return v * v.transpose();
	};
	
	const Mat cov() const {
		return scale * covariance / samples + epsilon * Mat::Identity();
	}
	
 public:
	AdaptiveMetropolisHastingsProposalDistribution() : 
		engine(rd()),
		normal_sample(std::normal_distribution<double>(0, 1)) {};
	
	void add(const Vec &sample) {
		// First sample.
		if (samples == 0) {
			sum = sample;
			samples++;
		} 
		
		// Second sample.
		else if (samples == 1) {
			Vec first = sum;
			sum += sample;
			samples++;
			covariance = outer(sample - sum / samples) + outer(first - sum / samples);
		}
		
		// Calculate covariance using a recurrence relation.
		else {
			Mat cov0 = samples * outer(sum / samples);
			sum += sample;
			samples++;
			Mat cov1 = samples * outer(sum / samples);
			covariance += cov0 - cov1 + outer(sample);
		}
		
		debug();
	};
	
	// Sample from a multivariate normal specified by the covariance.
	Vec next(const Vec &sample) const {
		if (samples < initial_period)
			return sample + multivariate_normal_sample(Mat::Identity());
		else
			return sample + multivariate_normal_sample(cov());
	};
	
	void debug() const {
		// For debug
		std::cerr << "Samples: " << samples << "; Variances: ";
		std::cerr << (cov().llt().matrixL() * Vec::Ones()).transpose() << std::endl;
	}
	
	// Convert functions to and from arrays.
	void add(const double *const a) {
		add(a2v(a));
	};
	void next(const double *const a, double *dest) {
		v2a(next(a2v(a)), dest);
	};
	
	void v2a(const Vec &v, double *dest) const {
		std::copy(v.data(), v.data() + dim, dest);
	};
	const Vec a2v(const double *const a) const {
		//return Eigen::Map<const Vec>(a);
		return Vec(a);
	};
};

#endif