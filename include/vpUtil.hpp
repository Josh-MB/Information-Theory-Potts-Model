#ifndef VP_TREE_UTIL_H
#define VP_TREE_UTIL_H
// A VP-Tree implementation, by Steve Hanov. (steve.hanov@gmail.com)
// http://stevehanov.ca/blog/index.php?id=130
// This implementation has been modified to use iterative approaches rather than recursive
// Additionally, we preallocate node storage rather than one at a time, since we intend to be
// using a large amount of nodes. Individual allocation results in a significant amount of overhead
// and reduced performance.
// Released to the Public Domain
// Based on "Data Structures and Algorithms for Nearest Neighbor Search" by Peter N. Yianilos

#include "defs.hpp"
#include "vpTree.hpp"

#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <queue>
#include <limits>
#include <stack>
#include <cmath>


template <unsigned int N>
struct VPoint {
	static const unsigned int length = N;
	double values[N];
	int idx;
};

/**
 * Distance functions for VPoints
 */
template <class T>
double VDistanceMaxnorm3D(const T &p1, const T &p2, const double L)
{
	double one = dwrap(p1.values[0] - p2.values[0], L);
	double two = dwrap(p1.values[1] - p2.values[1], L);
	double three = dwrap(p1.values[2] - p2.values[2], L);
	return std::max(one, std::max(two, three));
}

template <class T, unsigned int D1, unsigned int D2>
double VDistanceMaxnorm2D(const T &p1, const T &p2, const double L)
{
	double first = dwrap(p1.values[D1] - p2.values[D1], L);
	double second = dwrap(p1.values[D2] - p2.values[D2], L);
	return std::max(first, second);
}

template <class T, unsigned int D1, unsigned int D2>
double VDistanceEuclidean2D(const T &p1, const T &p2, const double L)
{
	double first = dwrap(p1.values[D1] - p2.values[D1], L);
	double second = dwrap(p1.values[D2] - p2.values[D2], L);
	return sqrt((first*first)+(second*second));
}

template <class T, unsigned int D1, unsigned int D2>
double VDistanceManhattan2D(const T &p1, const T &p2, const double L)
{
	double first = dwrap(p1.values[D1] - p2.values[D1], L);
	double second = dwrap(p1.values[D2] - p2.values[D2], L);
	return std::fabs(first) + std::fabs(second);
}

template <class T, unsigned int D1>
double VDistance1D(const T &p1, const T &p2, const double L)
{
	return dwrap(p1.values[D1] - p2.values[D1], L);
}

/**
 * used for calculating distance for (X,Y') dimensions in a (X,Y',W) space
 */
template <class T>
double VDistanceAllButLastDimension(const T &p1, const T &p2, const double L)
{
	double maxnorm = dwrap(p1.values[0] - p2.values[0], L);
	for(unsigned int i = 1; i < p1.length - 1; ++i) {
		maxnorm = std::max(maxnorm, dwrap(p1.values[i] - p2.values[i], L));
	}
	return maxnorm;
}

template <class T>
double VDistanceEuclideanAllButLastDimension(const T &p1, const T &p2, const double L)
{
	double dist = 0.0;
	for(unsigned int i = 0; i < p1.length - 1; ++i) {
		double d = dwrap(p1.values[i] - p2.values[i], L);
		dist += d*d;
	}
	return sqrt(dist);
}

template <class T>
double VDistanceManhattanAllButLastDimension(const T &p1, const T &p2, const double L)
{
	double dist = 0.0;
	for(unsigned int i = 0; i < p1.length - 1; ++i) {
		double d = dwrap(p1.values[i] - p2.values[i], L);
		dist += std::fabs(d);
	}
	return dist;
}

template <class T>
double VDistanceMaxnormAll(const T &p1, const T &p2, const double L)
{
	double maxnorm = dwrap(p1.values[0] - p2.values[0], L);
	for(unsigned int i = 1; i < p1.length; ++i) {
		maxnorm = std::max(maxnorm, dwrap(p1.values[i] - p2.values[i], L));
	}
	return maxnorm;
}

template <class T>
double VDistanceEuclideanAll(const T &p1, const T &p2, const double L)
{
	double dist = 0.0;
	for(unsigned int i = 0; i < p1.length; ++i) {
		double d = dwrap(p1.values[i] - p2.values[i], L);
		dist += d*d;
	}
	return sqrt(dist);
}

template <class T>
double VDistanceManhattanAll(const T &p1, const T &p2, const double L)
{
	double dist = 0.0;
	for(unsigned int i = 0; i < p1.length; ++i) {
		double d = dwrap(p1.values[i] - p2.values[i], L);
		dist += std::fabs(d);
	}
	return dist;
}

template <class T>
double VDistancePAll(const T &p1, const T &p2, const double L)
{
	double dist = 0.0;
	for(unsigned int i = 0; i < p1.length; ++i) {
		double d = dwrap(p1.values[i] - p2.values[i], L);
		dist += std::pow(std::fabs(d), 1.0/20.0);
	}
	return dist;
}

/**
 * Helpful type aliases
 */
using VpTreeEuclidean2D = VpTree< VPoint<2>, VDistanceEuclidean2D<VPoint<2>, 0, 1> >;
#endif //VP_TREE_UTIL_H