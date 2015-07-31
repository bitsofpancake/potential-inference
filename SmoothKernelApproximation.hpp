#ifndef SMOOTH_KERNEL_APPROXIMATION_HPP
#define SMOOTH_KERNEL_APPROXIMATION_HPP

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include "nanoflann.hpp"
#include "Particle.hpp"

using namespace boost::accumulators;

class SmoothKernelApproximation {
	std::vector<Particle> data;
	mutable accumulator_set<double, stats<tag::variance>> accumulator[2*dim];
	double bandwidth[2*dim];
	double vol = 1.0;
	
	double kernel(const double u) const;
	
 public:
	void add(const std::vector<Particle> others);
	void save();
	double operator()(const Particle a) const;
};



using namespace nanoflann;
class SmoothKernelApproximation2 {
		
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
	
	typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, KDAdaptor>, KDAdaptor, dim*2> KDTree;
	
	std::vector<Particle> data;
	mutable accumulator_set<double, stats<tag::variance>> accumulator[2*dim];
	double bandwidth[2*dim];
	double vol = 1.0;
	
	KDAdaptor adaptor;
	KDTree tree;
	
	double kernel(const double u) const;

 public:
	SmoothKernelApproximation2();
	void add(const std::vector<Particle> others);
	void save();
	double operator()(const Particle a) const;
};

#endif