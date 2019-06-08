#pragma once

#include "defs.hpp"
#include <vector>
#include <limits>
#include <algorithm>

class WangLandau
{
public:
	WangLandau(double const init_factor, double const factor_target, 
			   double const flatness_threshold, size_t const min_runs,
			   size_t const L, double const normalisation, size_t const normalisation_location);
	inline void visit(int E) { ++hist[E]; logG[E] += log_factor; }
	inline double tr_prob(int E1, int E2) const { return std::min(std::exp(logG[E1] - logG[E2]), 1.0); }
	bool update();
	bool isDone() const { return done; }
	size_t getIteration() const { return iteration; }
	std::vector<double> getLogG() const { return logG; }
	Hist getFinalHist() const { return final_hist; }
	Hist getHist() const { return hist; }
	//size_t getMinHistVal() const { return *std::min_element(hist.begin(), hist.end()); }
	double getFactor() const { return factor; }
	void setLogG(std::vector<double> const& newLogG, double const& newFactor);
	void setIteration(size_t it) { iteration = it; }
private:
	double factor;
	double log_factor;
	double factor_target;
	double flatness_threshold;
	double normalisation;
	size_t min_runs;
	size_t L;
	size_t normalisation_location;
	size_t iteration = 0;
	int sweepCount = 0;
	std::vector<double> logG;
	Hist hist;
	Hist final_hist;
	bool done = false;
	bool is_flat() const;
};