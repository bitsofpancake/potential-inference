#ifndef ADAPTIVE_METROPOLIS_HASTINGS_PROPOSAL_DISTRIBUTION_HPP
#define ADAPTIVE_METROPOLIS_HASTINGS_PROPOSAL_DISTRIBUTION_HPP

#include <Eigen/Dense>
using namespace Eigen;

// https://projecteuclid.org/download/pdf_1/euclid.bj/1080222083
template <int dim>
class AdaptiveMetropolisHastingsProposalDistribution {
	typedef Matrix<double, dim, dim> Mat;
	typedef Matrix<double, dim, 1> Vec;
	
	const double scale = 2.4 * 2.4 / dim; // Gelman, et al.
	const double epsilon = 0.01;
	const int initial_period = 10;
	
	Mat covariance;
	Vec sum;
	int samples = 0;
	
	std::random_device rd;
	std::mt19937 engine;
	std::normal_distribution<double> normal_sample(0, 1);
	
	inline Mat outer(const Vec v) const {
		return v * v.transpose();
	};
	
	Vec multivariate_normal_sample(const Mat cov) const {
		Vec random;
		for (int i = 0; i < dim; i++)
			random << normal_sample(engine);
		return sample + cov.llt() * random;
	}
	
 public:
	AdaptiveMetropolisHastingsProposalDistribution() : engine(rd()) {};
	
	void add(const Vec sample) {
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
			Mat cov0 = samples * outer(sum / samples)
			sum += sample;
			samples++;
			Mat cov1 = samples * outer(sum / samples);
			covariance += cov0 - cov1 + outer(sample);
		}
	};
	
	// Sample from a multivariate normal specified by the covariance.
	Vec next(const Vec sample) const {
		if (samples < initial_period)
			return multivariate_normal_sample(Mat::Identity());
		else
			return multivariate_normal_sample(scale * covariance / samples + epsilon * Mat::Identity());
	};
};

#endif