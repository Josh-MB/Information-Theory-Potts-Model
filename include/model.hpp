#pragma once

#include <random>
#include "defs.hpp"
#include "wanglandau.hpp"
#include "connectedSets.hpp"

void perfect_unordered_init(std::mt19937_64 &engine, std::vector<u8> &lattice, size_t const L);
void initialise(std::mt19937_64 &engine, std::vector<u8> &lattice, size_t const L, int const imode, double prop = 1.0);

constexpr int get_num_transitions();

// Initialise the transition probability look up table for specified temperature
std::vector<double> tp_table_init(const double T, char const potts_version);

// Performs an update on the lattice using glauber mechanics. Performs N site updates
// Note potts_version is unneeded here, due to the way calc_site_energy_delta and tpTable
// interact.
void update_glauber(std::mt19937_64 &engine,
					std::uniform_real_distribution<> &dist,
					std::vector<u8> &lattice, // Lattice data
					size_t const L,			// Lattice linear size
					size_t const N,			// Number of spin-flips per update
					std::vector<double> const& tpTable, // Transition probabilities
					int &current_energy
);

void update_swendsen_wang(std::mt19937_64& engine,
	std::uniform_real_distribution<>& dist,
	std::vector<u8>& lattice,
	size_t const L,
	size_t const N,
	double const T,
	int& current_energy,
	ConnectedSets& sets);

void update_glauber_and_record(std::mt19937_64 & engine, std::uniform_real_distribution<>& dist, std::vector<u8>& lattice, size_t const L, size_t const N, std::vector<double> const& tpTable, int &current_energy, FILE* fp);

void update_glauber_and_gte_hist(std::mt19937_64 & engine,
								 std::uniform_real_distribution<>& dist,
								 std::vector<u8>& lattice,
								 size_t const L,
								 size_t const N,
								 std::vector<double> const& tpTable,
								 int &current_energy,
								 std::vector<bool> const& petMask,
								 std::vector<Hist> & gte_hists,
								 Hist& gte_counts, 
								 std::vector<Hist> & mi_hists,
								 Hist& mi_counts,
								 std::vector<double> & magnetisation, 
								 Hist& mag_count, 
								 std::vector<Hist> & gte_binary_hists, 
								 Hist& gte_binary_counts, bool doScale = true);

void update_wang_landau(std::mt19937_64 &engine,
						std::uniform_real_distribution<> &dist,
						std::vector<u8> &lattice,
						WangLandau &wl,
						int &current_energy,
						size_t const L,
						size_t const N);

bool update_wang_landau_and_quantities(std::mt19937_64 & engine, std::uniform_real_distribution<>& dist, std::vector<u8>& lattice, WangLandau &wl, int & current_energy, size_t const L, size_t const N,
									   std::vector<Hist> & mi_hists, std::vector<double> &magnetisation, Hist& counts, int histCountLimit, std::vector<double> & inst_mi_vals, int const min_count_energy, double const sortThreshold,
									   Hist& interfacial_running_totals, Hist&interfacial_nums, std::vector<std::vector<u8>>& latticeQueue, int num_threads, std::vector<int> &energyQueue);


int calc_site_energy(std::vector<u8> const& lattice, size_t const x, size_t const y, size_t const L, char const potts_version);
int calc_site_energy_delta(std::vector<u8> const& lattice, u8 const newState, size_t const x, size_t const y, size_t const L);

int calculate_total_energy(std::vector<u8> const& lattice, size_t const L, char const potts_version);
double calc_magnetisation(std::vector<u8> const& lattice, size_t const L);
double calc_magnetisation_and_counts(std::vector<u8> const& lattice, size_t const L, Hist& count);
double critical_temp(char const potts_version);

int calc_action(std::vector<u8> const& lattice, size_t const L);
double calc_order(std::vector<u8> const & lattice, size_t const L);

double calculate_autocorrelation(std::vector<double> const& magnetisation, int dt, double mag_mean, double mag_var);
double calculate_relaxation_time(std::vector<double> const& magnetisation);

struct XY_site
{
	double cosTheta;
	double sinTheta;
	XY_site(double c, double s) : cosTheta(c), sinTheta(s) {};
private:
	friend double dot(const XY_site& lhs, const XY_site& rhs) {
		return lhs.cosTheta * rhs.cosTheta + lhs.sinTheta * rhs.sinTheta;
	}
};

// XY model update that uses Swendsen-Wang updating
double calc_energy_xy(std::vector<double>& lattice, std::vector<XY_site>& latticeVectors, std::vector<double>& siteEnergy, const size_t L);

double calc_site_energy_xy(std::vector<double>& lattice, std::vector<XY_site>& latticeVectors, XY_site const newState, size_t const x, size_t const y, const size_t L);

void update_glauber_xy(std::mt19937_64 & engine, std::uniform_real_distribution<>& dist,
			       std::vector<double>& swlattice,
				   std::vector<XY_site>& latticeVectors,
			       std::vector<double> &swSiteEnergy,
			       const double T, size_t const L, size_t const N, double &current_energy);

void update_swendsen_wang_xy(std::mt19937_64& engine, std::uniform_real_distribution<>& dist, std::vector<u8>& lattice, size_t const L, size_t const N, double const T, ConnectedSets& sets, const std::vector<double>& J);

void update_tomita(std::mt19937_64& engine, std::uniform_real_distribution<>& unit_dist, std::vector<double>& lattice, size_t const L, size_t const N, double const T, ConnectedSets& sets);
