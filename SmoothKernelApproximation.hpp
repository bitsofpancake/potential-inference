#ifndef SMOOTH_KERNEL_APPROXIMATION_HPP
#define SMOOTH_KERNEL_APPROXIMATION_HPP

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include "nanoflann.hpp"
#include "Particle.hpp"

class SmoothKernelApproximation_LinearSearch {
	typedef boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance>> VarianceAccumulator;

	std::vector<Particle> data;
	mutable VarianceAccumulator accumulator[2*dim];
	double bandwidth[2*dim];
	double vol = 1.0;
	
	double kernel(const double u) const;
	
 public:
	void add(const std::vector<Particle> others);
	void save();
	double operator()(const Particle a) const;
};

class SmoothKernelApproximation_KDTree {
		
	struct KDAdaptor {
		const std::vector<Particle> &data;
		const double (&bandwidth)[dim*2];
		KDAdaptor(std::vector<Particle> &data, double (&bandwidth)[dim*2]);
		inline size_t kdtree_get_point_count() const;
		inline double kdtree_distance(const double *b, const size_t j, size_t size) const;
		inline double kdtree_get_pt(const size_t i, int dim) const;
		
		template <class BoundingBox>
		bool kdtree_get_bbox(BoundingBox &x) const;
	};
	
	typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, KDAdaptor>, KDAdaptor, dim*2> KDTree;
	typedef boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance>> VarianceAccumulator;
	
	std::vector<Particle> data;
	mutable VarianceAccumulator accumulator[2*dim];
	double bandwidth[2*dim];
	double vol = 1.0;
	
	KDAdaptor adaptor;
	KDTree tree;
	
	double kernel(const double u) const;

 public:
	SmoothKernelApproximation_KDTree();
	void add(const std::vector<Particle> others);
	void save();
	double operator()(const Particle a) const;
};

typedef SmoothKernelApproximation_KDTree SmoothKernelApproximation;

#endif