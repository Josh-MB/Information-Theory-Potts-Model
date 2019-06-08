#include "../include/wanglandau.hpp"
#include "../include/utils.hpp"

WangLandau::WangLandau(double const init_factor, double const factor_target, double const flatness_threshold, size_t const min_runs, size_t const L, double const normalisation, size_t const normalisation_location) :
	factor(init_factor), factor_target(factor_target), flatness_threshold(flatness_threshold), normalisation(normalisation), min_runs(min_runs), L(L), normalisation_location(normalisation_location)
{
	log_factor = std::log(factor);
	logG.resize(L*L * 2+1, 1.);
	hist.resize(L*L * 2+1, 0u);
}

bool WangLandau::update()
{
	++sweepCount;
	//bool verbose = false;
	//if(sweepCount % 10000 == 0) {
		//printf("Sweeping... %d\n", sweepCount);
		//verbose = true;
	//}
	//printf("%d\n", sweepCount);

	if(sweepCount > min_runs && is_flat())
	{
		//printf("3\n");

		// Normalise DoS to get absolute density of states, or in other
		// words, g[groundState] = Q, where ground state is minimum energy (all same)
		double const offset = normalisation - logG[normalisation_location];
		std::transform(logG.begin(), logG.end(), logG.begin(), [offset](auto g) {return g + offset; });
		printf("Done iteration, factor = %f\n", factor);
		factor = std::sqrt(factor);
		//factor *= 0.5;
		log_factor = std::log(factor);
		final_hist = hist;
		if(factor < factor_target) {
			done = true;
		}
		clear_and_resize_vec(hist);
		sweepCount = 0;
		++iteration;
		return true;
	}
	return false;
}

void WangLandau::setLogG(std::vector<double> const & newLogG, double const& newFactor)
{
	if(newLogG.size() != logG.size()) PEEXIT("Trying to change size of logG when setting");
	factor = newFactor;
	log_factor = std::log(factor);
	logG = newLogG;
}

bool WangLandau::is_flat() const
{
	size_t hist_min = std::numeric_limits<size_t>::max();
	int num_bins = 0;
	size_t hist_sum = 0;
	
	size_t hist_size = hist.size();
	for(auto i = 0u; i < hist_size; ++i) 
	{
		if(hist[i] == 0) continue;
		++num_bins;
		hist_sum += hist[i];
		hist_min = std::min(hist[i], hist_min);
	}
	if(num_bins == 0) return false;
	double const mean_value = static_cast<double>(hist_sum) / num_bins;
	double const min_value_pct = hist_min / mean_value;

	//printf("Min val: %zu (%f), mean: %f\n", hist_min, min_value_pct, mean_value);
	return (min_value_pct >= flatness_threshold);
}
