#include "SmoothKernelApproximation.hpp"

using namespace boost::accumulators;
using namespace nanoflann;

// This implementation uses a simple linear search to compute the value of the approximation, in linear time.
double SmoothKernelApproximation_LinearSearch::kernel(const double u) const {
	if (u <= -1 || u >= 1)
		return 0;
	return 0.75 * (1 - u*u);
}

void SmoothKernelApproximation_LinearSearch::add(const std::vector<Particle> others) {
	data.insert(data.end(), others.begin(), others.end());
	for (const Particle &particle : others) {
		for (int i = 0; i < 2*dim; i++)
			accumulator[i](particle.coords[i]);
	}
};
	
void SmoothKernelApproximation_LinearSearch::save() {
	// Calculate bandwidth, via Scott's rule.
	for (int i = 0; i < 2*dim; i++) {
		bandwidth[i] = sqrt(variance(accumulator[i])) * pow(data.size(), -1.0 / (4 + 2 * dim));
		vol *= bandwidth[i];
	}
};
	
double SmoothKernelApproximation_LinearSearch::operator()(const Particle a) const {
	double value = 0.0;
	for (const Particle &b : data) {
		double u_sq = 0;
		for (int i = 0; i < 2*dim; i++)
			u_sq += pow((a.coords[i] - b.coords[i]) / bandwidth[i], 2);
		value += kernel(sqrt(u_sq));
		if (u_sq <= 1) {
		//	std::cout << sqrt(u_sq) << '\t' << b.coords[0] << '\t' << b.coords[1] << std::endl;
		}
	}
	return value / data.size() / vol;
};

// This implementation uses a k-d tree search algorithm to compute the value of the approximation, in logarithmic time.
double SmoothKernelApproximation_KDTree::kernel(const double u) const {
	if (u <= -1 || u >= 1)
		return 0;
	return 0.75 * (1 - u*u);
};
	
SmoothKernelApproximation_KDTree::SmoothKernelApproximation_KDTree() : 
	adaptor(data, bandwidth),
	tree(dim*2, adaptor) {};

void SmoothKernelApproximation_KDTree::add(const std::vector<Particle> &others) {
	points += others.size();
	data.reserve(points);
	
	for (const Particle &particle : others) {
		bool isnan = false;
		for (int i = 0; i < 2*dim; i++)
			if (!std::isfinite(particle.coords[i]))
				isnan = true;
		if (isnan)
			continue;
		
		data.push_back(particle);
		for (int i = 0; i < 2*dim; i++)
			accumulator[i](particle.coords[i]);
	}
};

void SmoothKernelApproximation_KDTree::save() {
	// Calculate bandwidth, via Scott's rule.
	for (int i = 0; i < 2*dim; i++) {
		bandwidth[i] = data.size() > 1 ? sqrt(variance(accumulator[i])) * pow(data.size(), -1.0 / (4 + 2 * dim)) : 1;
		vol *= bandwidth[i];
		
		// To-do: should we care if bandwidth[i] is infinite?
	}
	
	tree.buildIndex();
};

double SmoothKernelApproximation_KDTree::operator()(const Particle &a) const {	
	if (points == 0)
		return 0;

	Particle scaled;
	for (int i = 0; i < 2*dim; i++)
		scaled.coords[i] = a.coords[i] / bandwidth[i];
	
	SearchParams params;
	params.sorted = false;
	
	std::vector<std::pair<size_t, double>> neighbors;
	tree.radiusSearch(&(scaled.coords[0]), 1.0, neighbors, params);
	double value = 0.0;
	
	for (auto neighbor : neighbors) {
		//std::cout << sqrt(neighbor.second) << '\t' << data[neighbor.first].coords[0] << '\t' << data[neighbor.first].coords[1] << std::endl;
		value += kernel(sqrt(neighbor.second));
	}
		
	return value / points / vol;
};

SmoothKernelApproximation_KDTree::KDAdaptor::KDAdaptor(std::vector<Particle> &data, double (&bandwidth)[dim*2])
	: data(data), bandwidth(bandwidth) {};
inline size_t SmoothKernelApproximation_KDTree::KDAdaptor::kdtree_get_point_count() const {
	return data.size();
};
inline double SmoothKernelApproximation_KDTree::KDAdaptor::kdtree_distance(const double *b, const size_t j, size_t size) const {
	double u_sq = 0;
	for (int i = 0; i < 2*dim; i++)
		u_sq += pow(kdtree_get_pt(j, i) - b[i], 2);
	return u_sq;
};
inline double SmoothKernelApproximation_KDTree::KDAdaptor::kdtree_get_pt(const size_t j, int i) const {
	return data[j].coords[i] / bandwidth[i];
};
template <class BoundingBox>
bool SmoothKernelApproximation_KDTree::KDAdaptor::kdtree_get_bbox(BoundingBox &box) const {
	return false;
};