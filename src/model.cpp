#include "../include/model.hpp"
#include "../include/stats.hpp"
#include "../include/utils.hpp"
#include "../include/connectedSets.hpp"

#include <array>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <set>
#include <omp.h>

constexpr int num_transitions = 2 * neighbour_sites + 1;
constexpr int J = 1;

constexpr int get_num_transitions() {
	return num_transitions;
}

std::vector<double> tp_table_init(const double T, char const potts_version)
{
	std::vector<double> tpTable(num_transitions);
	double versionMultiplier = ((potts_version == 'b') + 1.0);
	for(int64_t n = -neighbour_sites; n <= neighbour_sites; ++n) {
		tpTable[n + neighbour_sites] = 1.0 / (1.0 + exp(versionMultiplier * n / T));
	}
	return tpTable;
}

u8 gen_state_exclude(std::mt19937_64 &engine, std::set<u8> exclude)
{
	if (exclude.size() > numStates) PEEXIT("Cannot exclude more values than there are states");
	int rangeLength = numStates - static_cast<int>(exclude.size());
	std::uniform_int_distribution<> rng(0, rangeLength - 1);
	u8 randomState = static_cast<u8>(rng(engine));
	for(auto e : exclude) {
		if(e > randomState) return randomState;
		++randomState;
	}
	return randomState;
}

void perfect_unordered_init(std::mt19937_64 &engine, std::vector<u8>& lattice, size_t const L)
{
	static_assert(numStates > 4, "Perfect unordered init only works for q>4");
	std::uniform_int_distribution<> rng(0, numStates - 1);

	// Make even chequerboard random
	for(auto y = 0u; y < L; ++y) {
		for(auto x = 0u; x < L; ++x) {
			auto sum = (x % 2) + (y % 2);
			if(sum % 2 == 0)
				lattice[y*L + x] = static_cast<u8>(rng(engine));
		}
	}
	// Fill in odd chequerboard such that no cell equals it neighbours
	std::set<u8> n;
	size_t const L1 = L - 1;
	for(auto y = 0u; y < L; ++y) {
		for(auto x = 0u; x < L; ++x) {
			auto sum = (x % 2) + (y % 2);
			if(sum % 2 == 1) {
				n.clear();
				n.insert(lattice[wrap_minus(y, L1)*L + x]);
				n.insert(lattice[wrap_plus(y, L1)*L + x]);
				n.insert(lattice[y*L + wrap_minus(x, L1)]);
				n.insert(lattice[y*L + wrap_plus(x, L1)]);
				lattice[y*L + x] = gen_state_exclude(engine, n);
			}
		}
	}
}

void initialise(std::mt19937_64 &engine, std::vector<u8>& lattice, size_t const L, int const imode, double prop)
{
	std::uniform_int_distribution<> rng(0, numStates-1);
	u8 staticState = 0;
	switch(imode)
	{
	case 0: // independently random
		for(auto &l : lattice) l = static_cast<u8>(rng(engine));
		break;
	case 1: // random but everything same
		staticState = static_cast<u8>(rng(engine));
		for(auto &l : lattice) l = staticState;
		break;
	case 2: // all state 0
		for(auto &l : lattice) l = staticState;
		break;
	case 4: // Init to P(E) valley, i.e <E>=-2.5
		// This isn't quite exact, but constructs a state with energy
		// equal to at least -2-2/L (and change from the random region)
		// This puts it near to <E>=-2.5, but favours the disordered side
		for(auto n = lattice.size(), i = decltype(n){0}; i < n; ++i) {
			if(i < (L* 4 * L / 8))
				lattice[i] = staticState;
			else
				lattice[i] = static_cast<u8>(rng(engine));
		}
		break;
	case 5: // Set prop of lattice to 0 and everything else random
	{
		std::uniform_real_distribution<> rng2(0.0, 1.0);
		for (auto& l : lattice) {
			if (rng2(engine) < prop) {
				l = staticState;
			}
			else {
				l = static_cast<u8>(rng(engine));
			}
		}
		break;
	}
	default:
		EEXIT("unknown initialisation mode");
	}
}

int calc_site_energy(std::vector<u8> const& lattice, size_t const x, size_t const y, size_t const L, char const potts_version)
{
	const size_t L1 = L - 1;
	size_t i = y*L + x;
	int delta_sum = (lattice[i] == lattice[wrap_minus(y, L1)*L + x]) +
		(lattice[i] == lattice[wrap_plus(y, L1)*L + x]) +
		(lattice[i] == lattice[y*L + wrap_minus(x, L1)]) +
		(lattice[i] == lattice[y*L + wrap_plus(x, L1)]);
	if(potts_version == 'a')
		return -J * (delta_sum);
	else
		return -J * (2 * delta_sum - 4);
}

// Calculates the delta in energy between two sites more efficiently than calculating each separately
// due to the decreased array look ups. Potts version is ignored since the delta is only used in indexing
// the transition table, as such, Version B results get halved immediately outside of this function and
// thus we skip the *2/2.
int calc_site_energy_delta(std::vector<u8> const& lattice, u8 const newState, size_t const x, size_t const y, size_t const L)
{
	size_t const L1 = L - 1;
	u8 const i = lattice[y*L + x];
	u8 const u = lattice[wrap_minus(y, L1)*L + x];
	u8 const d = lattice[wrap_plus(y, L1)*L + x];
	u8 const l = lattice[y*L + wrap_minus(x, L1)];
	u8 const r = lattice[y*L + wrap_plus(x, L1)];
	char delta_sum = ((i == u) + (i == d) + (i == r) + (i == l) - (newState == u) - (newState == d) - (newState == r) - (newState == l));
	return delta_sum;
}

void update_glauber(std::mt19937_64 & engine, std::uniform_real_distribution<>& dist, std::vector<u8>& lattice, size_t const L, size_t const N, std::vector<double> const& tpTable, int &current_energy)
{
	double const DL = static_cast<double>(L);
	std::uniform_int_distribution<> newStateDist(0, numStates-1);

	// N spin-flips per update using Glauber single-spin-flip updates
	// where P(S_i -> S_j) denotes the probability of random cell S flipping
	// to state j:
	// P(S_i -> S_j) = 1/(1+exp(S_3/T)),
	// where S_3 = S_2 - S_1, where S_2 is the site energy if the cell is updated to S_j
	// and S_1 is the site energy at the current state S_i
	for(size_t n = 0; n < N; ++n) {
		//Could speed this up by generating a single number
		//in the range 1:DL^2, and then calculate x,y from this
		//RNG took like 40% of calculation time from memory
		//so removing 1 call in the tight loop might speed
		//things up significantly
		const size_t x = static_cast<size_t>(dist(engine)*DL);
		const size_t y = static_cast<size_t>(dist(engine)*DL);
		
		const u8 newState = static_cast<u8>(newStateDist(engine));
		
		//Add neighbour sites to shift from [-n,n] (or [-2n,...-2,0,2,...,2n] for potts version b) range to [0,2n+1] range for indexing
		auto delta = calc_site_energy_delta(lattice, newState, x, y, L) + neighbour_sites;

		//Transition delta > U(0,1)
		if(tpTable[delta] > dist(engine)) {
			lattice[y*L + x] = newState;
			current_energy -= (delta - neighbour_sites);
		}
	}
}

void update_swendsen_wang(std::mt19937_64& engine, std::uniform_real_distribution<>& dist, std::vector<u8>& lattice, size_t const L, size_t const N, double const T, int& current_energy, ConnectedSets& sets)
{
	sets.clear();
	double const DL = static_cast<double>(L);
	size_t const L1 = L - 1;
	std::uniform_real_distribution<> fdist(0., 1.);
	double pFlip = 1 - exp(-2. / T);
	for (size_t y = 0; y < L; ++y) {
		for (size_t x = 0; x < L; ++x) {
			const std::array<size_t, 2> indices = { wrap_minus(y, L1) * L + x,
													y * L + wrap_minus(x, L1) };
			const size_t i_c = y * L + x;
			const u8 w_i = lattice[i_c];

			bool noAdj = true;
			for (size_t neighbour : indices)
			{
				// Not every particle in connected set will be part of virtual cluster
				if (w_i == lattice[neighbour] && fdist(engine) < pFlip)
				{
					sets.connect(i_c, neighbour);
					noAdj = false;
				}
			}
			// Didn't connect to anything so make sure it gets added to graph
			if (noAdj)
				sets.add(i_c);
		}
	}
	sets.squash();

	std::uniform_int_distribution<> newStateDist(0, numStates - 1);

	// Consider each virtual cluster once, and flip to random new state (potentially same as prev)
	for (auto& set : sets.sets) {
		const u8 newState = static_cast<u8>(newStateDist(engine));

		for (const int idx : set) {
			lattice[idx] = newState;
		}
	}
}

void update_glauber_and_record(std::mt19937_64 & engine, std::uniform_real_distribution<>& dist, std::vector<u8>& lattice, size_t const L, size_t const N, std::vector<double> const& tpTable, int &current_energy, FILE* fp)
{
	double const DL = static_cast<double>(L);
	size_t const L1 = L - 1;
	std::uniform_int_distribution<> newStateDist(0, numStates - 1);

	// N spin-flips per update using Glauber single-spin-flip updates
	// where P(S_i -> S_j) denotes the probability of random cell S flipping
	// to state j:
	// P(S_i -> S_j) = 1/(1+exp(S_3/T)),
	// where S_3 = S_2 - S_1, where S_2 is the site energy if the cell is updated to S_j
	// and S_1 is the site energy at the current state S_i
	size_t numSame = 0;
	for(size_t n = 0; n < N; ++n) {
		//Could speed this up by generating a single number
		//in the range 1:DL^2, and then calculate x,y from this
		//RNG took like 40% of calculation time from memory
		//so removing 1 call in the tight loop might speed
		//things up significantly
		const size_t x = static_cast<size_t>(dist(engine)*DL);
		const size_t y = static_cast<size_t>(dist(engine)*DL);

		const u8 newState = static_cast<u8>(newStateDist(engine));

		//Add neighbour sites to shift from [-n,n] (or [-2n,...-2,0,2,...,2n] for potts version b) range to [0,2n+1] range for indexing
		auto delta = calc_site_energy_delta(lattice, newState, x, y, L) + neighbour_sites;

		//Transition delta > U(0,1)
		u8 x_i = lattice[y*L + x];// lattice[y*L + x];

		if(tpTable[delta] > dist(engine)) {
			lattice[y*L + x] = newState;
			current_energy -= (delta - neighbour_sites);
		}

		u8 w_i = lattice[y*L + x];
		numSame += (x_i == w_i);
		u8 s_u = lattice[wrap_minus(y, L1)*L + x]; //This is new & old lattice (only x/w changes)
		u8 s_d = lattice[wrap_plus(y, L1)*L + x];
		u8 s_l = lattice[y*L + wrap_minus(x, L1)];
		u8 s_r = lattice[y*L + wrap_plus(x, L1)];

		BIN_WRITE(w_i, u8, 1, fp);
		BIN_WRITE(x_i, u8, 1, fp);
		BIN_WRITE(s_r, u8, 1, fp);
		BIN_WRITE(s_l, u8, 1, fp);
		BIN_WRITE(s_d, u8, 1, fp);
		BIN_WRITE(s_u, u8, 1, fp);
	}
	printf("Num same (per flip): %zu\n", numSame);
}

void update_glauber_and_gte_hist(std::mt19937_64 & engine, std::uniform_real_distribution<>& dist, std::vector<u8>& lattice, size_t const L, size_t const N, std::vector<double> const& tpTable, int &current_energy, std::vector<bool> const& petMask, std::vector<Hist> & gte_hists, Hist& gte_counts, std::vector<Hist> & mi_hists, std::vector<size_t>& mi_counts, std::vector<double> & magnetisation, Hist& mag_count, std::vector<Hist> & gte_binary_hists, Hist& gte_binary_counts, bool doScale)
{
	double const DL = static_cast<double>(L);
	size_t const L1 = L - 1;
	std::uniform_int_distribution<> newStateDist(0, numStates - 1);
	int const len = numStates * numStates * numSiteEnergy;
	
	// Calculate current magnetisation array (Just a histogram of all lattice
	// states).
	Hist count(numStates, 0);
	for(auto l : lattice) {
		++count[l];
	}

	// N spin-flips per update using Glauber single-spin-flip updates
	// where P(S_i -> S_j) denotes the probability of random cell S flipping
	// to state j:
	// P(S_i -> S_j) = 1/(1+exp(S_3/T)),
	// where S_3 = S_2 - S_1, where S_2 is the site energy if the cell is updated to S_j
	// and S_1 is the site energy at the current state S_i
	//int f = 0, p = 0;
	for(size_t n = 0; n < N; ++n) {
		//Could speed this up by generating a single number
		//in the range 1:DL^2, and then calculate x,y from this
		//RNG took like 40% of calculation time from memory
		//so removing 1 call in the tight loop might speed
		//things up significantly
		const size_t x = static_cast<size_t>(dist(engine)*DL);
		const size_t y = static_cast<size_t>(dist(engine)*DL);
		auto idx = y * L + x;

		const u8 newState = static_cast<u8>(newStateDist(engine));

		//Add neighbour sites to shift from [-n,n] (or [-2n,...-2,0,2,...,2n] for potts version b) range to [0,2n+1] range for indexing
		auto delta = calc_site_energy_delta(lattice, newState, x, y, L) + neighbour_sites;

		u8 x_i = lattice[idx];
		auto old_energy = current_energy;
		//Transition delta > U(0,1)
		if(tpTable[delta] > dist(engine)) {
			lattice[idx] = newState;
			current_energy -= (delta - neighbour_sites);
		}
		u8 w = lattice[idx];

		// Update magnetisation histogram
		--count[x_i];
		++count[w];

		//++f;
		//auto scaledEnergy = old_energy;// scaleEnergyValue(old_energy, L);
		auto scaledEnergy = doScale ? scaleEnergyValue(old_energy, L) : old_energy;
		if(scaledEnergy > petMask.size())
			PEEXIT("out of bounds: %d\n", old_energy);
		if(petMask[scaledEnergy]) {
			u8 s_i = x_i;// lattice[y*L + x];
			u8 s_u = lattice[wrap_minus(y, L1)*L + x]; //This is new & old lattice (only x/w changes)
			u8 s_d = lattice[wrap_plus(y, L1)*L + x];
			u8 s_l = lattice[y*L + wrap_minus(x, L1)];
			u8 s_r = lattice[y*L + wrap_plus(x, L1)];
			int e_x = (s_i == s_u) + (s_i == s_d) + (s_i == s_l) + (s_i == s_r);
			//int e_w = (w == s_u) + (w == s_d) + (w == s_l) + (w == s_r);
			int bin = e_x;// delta; // e_x - e_w + 4;
			int hist_idx = w + numStates * (x_i + numStates * bin);
			if(hist_idx >= len || hist_idx < 0) PEEXIT("Bad idx: %d, delta: %d\n", hist_idx, delta);// PEEXIT("2: %d, bin: %d, e_x %d, e_w, %d", idx, bin, e_x, e_w);
			++gte_hists[scaledEnergy][hist_idx];
			++gte_counts[scaledEnergy];

			int binary_hist_idx = w + numStates * (x_i +
													numStates * ((x_i == s_u) +
													2 * ((x_i == s_d) +
														 2 * ((x_i == s_l) +
															  2 * ((x_i == s_r))))));
			++gte_binary_hists[scaledEnergy][binary_hist_idx];
			++gte_binary_counts[scaledEnergy];

			//++p;

			/*++mi_hists[scaledEnergy][s_i + numStates * s_r];
			++mi_hists[scaledEnergy][s_i + numStates * s_l];
			++mi_hists[scaledEnergy][s_i + numStates * s_d];
			++mi_hists[scaledEnergy][s_i + numStates * s_u];
			mi_counts[scaledEnergy] += 4;*/

			// Work out most plentiful state, and calc magnetisation using it
			// where M = [(q*N_max/N)-1]/(q-1) according to binder81
			auto max = static_cast<double>(*std::max_element(count.begin(), count.end()));
			auto toAdd = ((numStates * max / (L*L)) - 1.0) / (numStates - 1.0);
			if(toAdd + magnetisation[scaledEnergy] < magnetisation[scaledEnergy]) {
				PEEXIT("magnetisation overflow: %f, %d, %f\n", magnetisation[scaledEnergy], scaledEnergy, toAdd);
			}
			magnetisation[scaledEnergy] += toAdd;
			++mag_count[scaledEnergy];
		}
	}
	//printf("One sweep: %zu, flipped: %d, pet: %d\n", N, f, p);
}

void update_wang_landau(std::mt19937_64 & engine, std::uniform_real_distribution<>& dist, std::vector<u8>& lattice, WangLandau &wl, int & current_energy, size_t const L, size_t const N)
{
	double const DL = static_cast<double>(L);
	std::uniform_int_distribution<> newStateDist(0, numStates - 1);

	for(size_t n = 0; n < N; ++n) {
		const size_t x = static_cast<size_t>(dist(engine)*DL);
		const size_t y = static_cast<size_t>(dist(engine)*DL);

		const u8 newState = static_cast<u8>(newStateDist(engine));

		// Negated as calc_site_energy was developed for glauber update
		// Doubled to account for the fact that neighbours will also update energy
		auto delta = -calc_site_energy_delta(lattice, newState, x, y, L);

		auto prob = wl.tr_prob(current_energy, current_energy + delta);
		if(prob > dist(engine)) {
			lattice[y*L + x] = newState;
			//Update the total energy
			current_energy += delta;
		}
		wl.visit(current_energy);
	}
}

bool update_wang_landau_and_quantities(std::mt19937_64 & engine, std::uniform_real_distribution<>& dist, std::vector<u8>& lattice, WangLandau &wl, int & current_energy, size_t const L, size_t const N,
									   std::vector<Hist> & mi_hists, std::vector<double> &magnetisation, Hist& counts, int histCountLimit, std::vector<double> & inst_mi_vals,
									   int const min_count_energy, double const sortThreshold, Hist&interfacial_running_totals, Hist&interfacial_nums, std::vector<std::vector<u8>>& latticeQueue, int num_threads, std::vector<int> &energyQueue)
{
	double const DL = static_cast<double>(L);
	bool retVal = false;
	std::uniform_int_distribution<> newStateDist(0, numStates - 1);
	Hist inst_mi(numStates*numStates, 0);
	int threadCounter = 0;
	for(size_t n = 0; n < N; ++n) {
		const size_t x = static_cast<size_t>(dist(engine)*DL);
		const size_t y = static_cast<size_t>(dist(engine)*DL);

		const u8 newState = static_cast<u8>(newStateDist(engine));

		// Negated as calc_site_energy was developed for glauber update
		// Doubled to account for the fact that neighbours will also update energy
		auto delta = -calc_site_energy_delta(lattice, newState, x, y, L);

		auto prob = wl.tr_prob(current_energy, current_energy + delta);
		if(prob > dist(engine)) {
			lattice[y*L + x] = newState;
			//Update the total energy
			current_energy += delta;
		}
		wl.visit(current_energy);

		if(counts[current_energy] < histCountLimit) {
			Hist freq(numStates, 0);
			// Accumulate magnetisation values. Can use counts for averaging
			//reenable after interface testing
			double mag = calc_magnetisation_and_counts(lattice, L, freq);
			magnetisation[current_energy] += mag;

			latticeQueue[threadCounter] = lattice;
			energyQueue[threadCounter] = current_energy;
			++threadCounter;

			Hist length_sum(num_threads);
			Hist length_sizes(num_threads);

			if(threadCounter == num_threads) {
				//printf("Thread counter hit\n");
#pragma omp parallel for
				for(int i = 0; i < num_threads; ++i) {
					auto lengths = calcInterfaceAllLengths(latticeQueue[i], static_cast<int>(L));
					length_sum[i] = std::accumulate(lengths.begin(), lengths.end(), size_t{}, std::plus<>());
					length_sizes[i] += lengths.size();
				}
				//printf("Accumulated\n");
				for(int i = 0; i < num_threads; ++i) {
					current_energy = energyQueue[i];

					interfacial_running_totals[current_energy] += length_sum[i];
					interfacial_nums[current_energy] += length_sizes[i];
				}
				threadCounter = 0;
			}
			//auto lengths = calcInterfaceAllLengths(lattice, L);
			//interfacial_running_totals[current_energy] += std::accumulate(lengths.begin(), lengths.end(), 0u, std::plus<>());
			//interfacial_nums[current_energy] += lengths.size();

			//double mean_length = std::accumulate(lengths.begin(), lengths.end(), 0., std::plus<double>()) / lengths.size();
			//interfacial_lengths[current_energy]->insert(interfacial_lengths[current_energy]->end(), lengths.begin(), lengths.end());

			// reenable after interface testing
			Hist sort_order(numStates, 0);
			std::iota(sort_order.begin(), sort_order.end(), 0);

			if(mag > sortThreshold)
			{
				auto temp_sort_order = sort_indexes(freq);
				for(int i = 0; i < numStates; ++i) {
					sort_order[temp_sort_order[i]] = i;
				}

			}
			// Accumulate MI histograms
			for(auto y1 = 0u; y1 < L; ++y1) {
				for(auto x1 = 0u; x1 < L; ++x1) {
					int64_t s_i = sort_order[lattice[y1*L + x1]];
					int64_t s_r = sort_order[lattice[y1*L + wrap_plus(x1, L - 1)]];
					int64_t s_d = sort_order[lattice[wrap_plus(y1, L - 1)*L + x1]];
					int64_t s_l = sort_order[lattice[y1*L + wrap_minus(x1, L - 1)]];
					int64_t s_u = sort_order[lattice[wrap_minus(y1, L - 1)*L + x1]];

					++mi_hists[current_energy][s_i + numStates * s_r];
					++mi_hists[current_energy][s_i + numStates * s_l];
					++mi_hists[current_energy][s_i + numStates * s_d];
					++mi_hists[current_energy][s_i + numStates * s_u];

					//++inst_mi[s_i + numStates * s_r];
					//++inst_mi[s_i + numStates * s_l];
					//++inst_mi[s_i + numStates * s_d];
					//++inst_mi[s_i + numStates * s_u];
				}
			}
			++counts[current_energy];

			if(current_energy == min_count_energy)
				retVal = true;
		}
	}

	return retVal;
}

int calculate_total_energy(std::vector<u8> const & lattice, size_t const L, char const potts_version)
{
	int energy = 0;
	for(auto y = 0u; y < L; ++y) {
		for(auto x = 0u; x < L; ++x) {
			energy += calc_site_energy(lattice, x, y, L, potts_version);
		}
	}
	return energy;
}

double calc_magnetisation(std::vector<u8> const& lattice, size_t const L)
{
	Hist count(numStates, 0);
	return calc_magnetisation_and_counts(lattice, L, count);
}

double calc_magnetisation_and_counts(std::vector<u8> const& lattice, size_t const L, Hist& count)
{
	for(auto l : lattice) {
		++count[l];
	}
	auto max = static_cast<double>(*std::max_element(count.begin(), count.end()));
	return ((numStates * max / (L*L)) - 1.0) / (numStates - 1.0);
}

double critical_temp(char const potts_version)
{
	// Double critical temperature if it is version b
	return ((potts_version == 'b') + 1.0) / std::log(1 + std::sqrt(numStates));
}

int calc_action(std::vector<u8> const & lattice, size_t const L)
{
	auto action = 0;
	size_t const L1 = L - 1;
	for(auto y = 0u; y < L; ++y) {
		for(auto x = 0u; x < L; ++x) {
			u8 s_i = lattice[y*L + x];
			u8 s_r = lattice[y*L + wrap_plus(x, L1)];
			u8 s_d = lattice[wrap_plus(y, L1)*L + x];
			action += (s_i == s_r) + (s_i == s_d);
		}
	}
	return action;
}

double calc_order(std::vector<u8> const & lattice, size_t const L)
{
	auto order = 0.;
	for(auto y = 0u; y < L; ++y) {
		for(auto x = 0u; x < L; ++x) {
			order += calc_site_energy(lattice, x, y, L, 'a') / 4.0;
		}
	}
	return order;
}

double calculate_autocorrelation(std::vector<double> const & magnetisation, int dt, double mag_mean, double mag_var)
{
	int64_t T = static_cast<int64_t>(magnetisation.size());
	double sum = 0;
	for(int64_t i = 0; i < T-dt; ++i) {
		sum += (magnetisation[i + dt] - mag_mean) * (magnetisation[i] - mag_mean);
	}
	return (sum / (T - dt)) / mag_var;
}

double calculate_relaxation_time(std::vector<double> const & magnetisation)
{
	double mag_mean = std::accumulate(magnetisation.begin(), magnetisation.end(), 0.) / magnetisation.size();
	double mag_var = std::accumulate(magnetisation.begin(), magnetisation.end(), 0.,
									 [mag_mean](double total, double x) {return total + (x - mag_mean)*(x - mag_mean); }) / magnetisation.size();
	double tau = 0.;
	for(int i = 0; i < 100; ++i) {
		double cfunc = calculate_autocorrelation(magnetisation, i, mag_mean, mag_var);
		if(cfunc < 0) break;
		tau += cfunc;
	}
	return tau;
}
