#include <vector>
#include <cstddef>
#include <random>
#include "defs.hpp"

/**
 * Calculates MI histogram from given lattice. Assumes hist is pre-zeroed.
 * Counts are half of actual. 
 */
void MI_histogram(std::vector<u8> const& lattice, size_t const L, Hist& hist);

/**
* Calculates TE histogram from given lattice. Assumes hist is pre-zeroed.
*/
void TE_histogram(std::vector<u8> const& lattice, std::vector<u8> const& lattice_old, size_t const L, Hist& hist);

/**
* Calculates GTE histogram from given lattice. Assumes hist is pre-zeroed.
*/
void GTE_histogram(std::vector<u8> const& lattice, std::vector<u8> const& lattice_old, size_t const L, Hist& hist);
void GTE_reduced_histogram(std::vector<u8> const& lattice, std::vector<u8> const& lattice_old, size_t const L, Hist& hist);

void all_histograms(std::vector<u8> const& lattice, std::vector<u8> const& lattice_old, size_t const L, Hist& mi, Hist& te, Hist& gte, Hist&binary_gte, Hist& reduced_gte, bool sortStates, std::mt19937_64 & engine);

double calc_MI(Hist const& hist);
double calc_GTE_from_hist(Hist const& hist);
double calc_GTE_reduced_from_hist(Hist const & hist);
double calc_GTE_binary_from_hist(Hist const & hist);
double entropy(std::vector<double> const& dist);
void normalise_hist(Hist const& hist, std::vector<double> & normalised_hist);

double kl_divergence(Hist const& baseHist, Hist const& targetHist);

/**
 * Scales energy value at current L to the corresponding energy level at L=32
 */
int scaleEnergyValue(int const energy, size_t const L);

std::vector<int> calcInterfaceAllLengths(std::vector<u8> const& lattice, int const L);
std::vector<int> calcInterfaceAllLengthsAlt(std::vector<u8> const& lattice, int const L);
