#include "SmoothKernelApproximation.hpp"

double SmoothKernelApproximation::kernel(const double u) const {
	if (u <= -1 || u >= 1)
		return 0;
	return 0.75 * (1 - u*u);
}

void SmoothKernelApproximation::add(const std::vector<Particle> others) {
	data.insert(data.end(), others.begin(), others.end());
	for (const Particle &particle : others) {
		for (int i = 0; i < 2*dim; i++)
			accumulator[i](particle.coords[i]);
	}
};
	
void SmoothKernelApproximation::save() {
	// Calculate bandwidth, via Scott's rule.
	for (int i = 0; i < 2*dim; i++) {
		bandwidth[i] = sqrt(variance(accumulator[i])) * pow(data.size(), -1.0 / (4 + 2 * dim));
		vol *= bandwidth[i];
	}
};
	
double SmoothKernelApproximation::operator()(const Particle a) const {
	double value = 0.0;			
	for (const Particle &b : data) {
		double u_sq = 0;
		for (int i = 0; i < 2*dim; i++)
			u_sq += pow((a.coords[i] - b.coords[i]) / bandwidth[i], 2);
		value += kernel(sqrt(u_sq));
	}
	return value / data.size() / vol;
};


SmoothKernelApproximation2::KDAdaptor::KDAdaptor(std::vector<Particle> &data, double (&bandwidth)[dim*2])
	: data(data), bandwidth(bandwidth) {};
inline size_t SmoothKernelApproximation2::KDAdaptor::kdtree_get_point_count() const {
	return data.size();
};
inline double SmoothKernelApproximation2::KDAdaptor::kdtree_distance(const double *b, const size_t j, size_t size) const {
	double u_sq = 0;
	for (int i = 0; i < 2*dim; i++)
		u_sq += pow((data[j].coords[i] - b[i]) / bandwidth[i], 2);
	return u_sq;
};
inline double SmoothKernelApproximation2::KDAdaptor::kdtree_get_pt(const size_t i, int dim) const {
	return data[i].coords[dim];
};
template <class BoundingBox>
bool SmoothKernelApproximation2::KDAdaptor::kdtree_get_bbox(BoundingBox &x) const {
	return false;
};

using namespace nanoflann;
double SmoothKernelApproximation2::kernel(const double u) const {
	if (u <= -1 || u >= 1)
		return 0;
	return 0.75 * (1 - u*u);
};
	
SmoothKernelApproximation2::SmoothKernelApproximation2() : 
	adaptor(data, bandwidth),
	tree(dim*2, adaptor) {};

void SmoothKernelApproximation2::add(const std::vector<Particle> others) {
	data.insert(data.end(), others.begin(), others.end());
	for (const Particle &particle : others) {
		for (int i = 0; i < 2*dim; i++)
			accumulator[i](particle.coords[i]);
	}
};

void SmoothKernelApproximation2::save() {
	// Calculate bandwidth, via Scott's rule.
	for (int i = 0; i < 2*dim; i++) {
		bandwidth[i] = sqrt(variance(accumulator[i])) * pow(data.size(), -1.0 / (4 + 2 * dim));
		vol *= bandwidth[i];
	}
	
	tree.buildIndex();
};

double SmoothKernelApproximation2::operator()(const Particle a) const {
	SearchParams params;
	std::vector<std::pair<size_t, double>> neighbors;
	
	tree.radiusSearch(&(a.coords[0]), 1.0, neighbors, params);
	double value = 0.0;
	for (auto neighbor : neighbors)
		value += kernel(sqrt(neighbor.second));
		
	return value / data.size() / vol;
};