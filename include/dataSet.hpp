#ifndef DATASET_HPP
#define DATASET_HPP

#include "vpUtil.hpp"
#include "psi.hpp"
#include "defs.hpp"
#include "vpTree.hpp"
#include "utils.hpp"

#include <omp.h>
#include <random>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fmt/format.h>

/**
 * Holds data and functions for calculating dimensional quantites
 * needed by nearest-neighbour estimators (e.g., hyper-ball radii
 * and neighbour counting).
 * Templated on the number of dimensions stored in the dataset
 */
template <unsigned int D>
class DataSet
{
public:
	// Sets up data set holding N elements, using k nearest neighbours
	// for determining hyper-ball sizes
	DataSet(int const N, int const k, std::mt19937_64* const rng = nullptr) :
		N(N), k(k), rng(rng)
	{
		ballsizes.resize(N);
		points.resize(N);
	}

	// Get the elements being stored
	VPoint<D> const* getPoints() const {
		return points.data();
	}
	// Get the calculated hyper-ball sizes
	double const* getBallsizes() const {
		return ballsizes.data();
	}

	// Fill dimension d with data stream
	void setData(double const* const data, unsigned int const d) {
		if(d >= 0 && d < dim)
		{
			for(int i = 0; i < N; ++i) {
				points[i].values[d] = data[i];
				points[i].idx = i;
			}
		}
	}
	// Set element i of dimension d with data
	inline void setData(double const data, int i, unsigned int d) {
		points[i].values[d] = data;
		points[i].idx = i;
	}

	inline int getN() const { return N; }
	inline int getK() const { return k; }

	// Calculate hyper-ball radii of dataset, using max-norm
	// distances
	void calcBallsizes() {
		VpTree<VPoint<D>, VDistanceMaxnormAll<VPoint<D> > > tree(TWOPI, rng);
		tree.create(points);

#pragma omp parallel for
		for(int i = 0; i < N; ++i) {
			VPoint<D> tmp;
			tree.find_kth_neighbour(points[i], (int)(k + 1), tmp, ballsizes[points[i].idx]);
		}
	}

	// Calculate the average ballsize of dataset
	double averageBallsize() {
		double total = std::accumulate(ballsizes.begin(), ballsizes.end(), 0.0,
			[](double const& acc, double const& x) {return acc + log(x * 2); });
		return (total/N);
	}

	// Calculate <psi(N_d(i)> for averaged every point along dimension D1
	// using the precalculated ballsizes
	template <unsigned int D1>
	double countNeighboursOneDim() {
		int d = D1;
		double total = 0.0;

		std::vector<double> sortedMarginal(N);
		for(int i = 0; i < N; ++i)
			sortedMarginal[i] = points[i].values[d];
		std::sort(sortedMarginal.begin(), sortedMarginal.end());
		int numThreads = 1;
#pragma omp parallel
		numThreads = omp_get_num_threads();
		std::vector<double> totals(numThreads);
#pragma omp parallel for
		for(int i = 0; i < N; ++i) {
			int const threadNum = omp_get_thread_num();
			totals[threadNum] += binarySearchRange(sortedMarginal, points[i].values[d], ballsizes[points[i].idx]);
		}
		
		for(int i = 0; i < numThreads; ++i) {
			total += totals[i];
		}
		return (total / N);
	}

	// Calculate <psi(N_d(i)> for every point along dimensions D1, D2
	// using the precalculated ballsizes
	template <unsigned int D1, unsigned int D2>
	double countNeighboursTwoDim() {
		VpTree<VPoint<D>, VDistanceMaxnorm2D<VPoint<D>, D1, D2> > tree;
		tree.create(points);
		double total = 0.0;
		int numThreads = 1;
#pragma omp parallel
		numThreads = omp_get_num_threads();
		std::vector<double> totals(numThreads);
#pragma omp parallel for
		for(int i = 0; i < N; ++i) {
			int const threadNum = omp_get_thread_num();
			int const n = tree.fr_search(points[i], ballsizes[points[i].idx]);
			totals[threadNum] += digamma(n);
		}

		for(int i = 0; i < numThreads; ++i) {
			total += totals[i];
		}

		return (total / N);
	}

	// Calculate <psi(N_d(i)> for every point along all dimensions except
	// the last using the precalculated ballsizes.
	// This is used for calculating N_xy, where y is multidimensional,
	// assuming dataset stores [X, Y', W]
	double countNeighboursNDim() {
		VpTree<VPoint<D>, VDistanceAllButLastDimension<VPoint<D> > > tree;
		tree.create(points);
		double total = 0.0;
		int numThreads = 1;
#pragma omp parallel
		numThreads = omp_get_num_threads();
		std::vector<double> totals(numThreads);
#pragma omp parallel for
		for(int i = 0; i < N; ++i) {
			int const threadNum = omp_get_thread_num();
			int const n = tree.fr_search(points[i], ballsizes[points[i].idx]);
			totals[threadNum] += digamma(n);
		}
		for(int i = 0; i < numThreads; ++i) {
			total += totals[i];
		}
		return (total / N);
	}

	static unsigned int const dim = D;
private:
	std::vector<VPoint<D> > points;
	std::vector<double> ballsizes;
	int const N;
	int const k;
	std::mt19937_64* const rng;

	// Perform binary search to find the number of elements, n, in sortedMarginal
	// that exist in the range base +- ballsize, and returns psi(n)
	double binarySearchRange(std::vector <double> const& sortedMarginal, double const base, double const ballsize) const
	{
		int baseIdx = binarySearch(sortedMarginal, base);

		double nextballsize = nextafter(ballsize, ballsize - 1);

		double lowerSearch = base - nextballsize;
		lowerSearch = rwrap(lowerSearch, TWOPI);

		int lowerIdx = binarySearch(sortedMarginal, lowerSearch);
		if (lowerIdx == -1) PEEXIT("LowerIdx = -1");

		double upperSearch = base + nextballsize;
		upperSearch = rwrap(upperSearch, TWOPI);

		int upperIdx = binarySearch(sortedMarginal, upperSearch);
		if (upperIdx == -1) PEEXIT("UpperIdx = -1");

		//Move lowerIdx to first element that succeeds (base - lower < ballsize)
		if (dwrap(sortedMarginal[lowerIdx] - base, TWOPI) < ballsize) {
			while (lowerIdx >= 0 && lowerIdx != upperIdx) {
				int nextIdx = (lowerIdx - 1 + N) % N;
				if (dwrap(sortedMarginal[nextIdx] - base, TWOPI) < ballsize) {
					lowerIdx = nextIdx;
				}
				else
					break;
			}
		}
		else {
			while (lowerIdx < N && lowerIdx != upperIdx) {
				lowerIdx = (lowerIdx + 1) % N;
				if (dwrap(sortedMarginal[lowerIdx] - base, TWOPI) < ballsize)
					break;
			}
		}

		//Move upperIdx to first element that fails (upper - base < ballsize)
		if (dwrap(sortedMarginal[upperIdx] - base, TWOPI) < ballsize) {
			while (upperIdx < N && lowerIdx != upperIdx) {
				upperIdx = (upperIdx + 1) % N;
				if (dwrap(sortedMarginal[upperIdx] - base, TWOPI) >= ballsize)
					break;
			}
		}
		else {
			while (upperIdx >= 0 && lowerIdx != upperIdx) {
				int nextIdx = (upperIdx - 1 + N) % N;
				if (dwrap(sortedMarginal[nextIdx] - base, TWOPI) < ballsize)
					break;
				else {
					upperIdx = nextIdx;
				}
			}
		}

		int n = upperIdx - lowerIdx;

		if (n < 0)
			n += N;

		//It is not possible for n = 0, since we count the point itself as well (instead of just adding 1 like kraskov. Same thing)
		//It is possible for n = numPts, in the case when there is large spread in one dimension, but not in the other. The ballsize
		//will be based off the large spread dimension and may engulf the entire spread in the lower spread dimension. Rare but possible.
		if (n == 0)
		{
			//If element next to base is within the bounds, but n == 0, then all elements must be within the bounds
			if (upperIdx == lowerIdx && dwrap(sortedMarginal[lowerIdx] - base, TWOPI) == ballsize) {
				n = N - 1;
			}
			else if (dwrap(sortedMarginal[(baseIdx + 1) % N] - base, TWOPI) < ballsize) {
				n = N;

			}
			else
				n = 1;
		}
		return digamma(n);
	}

	/**
	* \brief Standard binary search
	* If the query point does not exist, will return the index at which that point would be inserted. Allows bounds to be found without needing points with those values to exist
	* \param data - Sorted array to search
	* \param query - Value searching for
	* \return Lowest index occurence of query in data, or if does not exist, index where it would be inserted.
	*/
	int binarySearch(std::vector<double> const& data, double const query) const {
		int imin = 0;
		int imax = (int)data.size() - 1; //Inclusive indices

		while(imin < imax)
		{
			int imid = (imin + imax) >> 1;

			if(imid >= imax) return -1;

			if(data[imid] < query)
				imin = imid + 1;
			else
				imax = imid;
		}

		if(imax == imin)// && (data[imin] == query))
			return imin;
		return -1; //Need to return position would insert into
	}
};

#endif //DATASET_HPP
