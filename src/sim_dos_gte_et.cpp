#include "../include/model.hpp"
#include "../include/defs.hpp"
#include "../include/stats.hpp"
#include "../include/wanglandau.hpp"
#include "../include/version.hpp"
#include "../include/utils.hpp"
#include "../include/connectedSets.hpp"

#include <clara.hpp>
#include <fmt/format.h>
#include <vector>
#include <algorithm>
#include <random>
#include <omp.h>
#include <bitset>


/**
* Used for calculating a profile of GTE(E,T) for a given temperature over
* all energy values, where range of E is determined by lattice size.
* This can then be used with the density of states to calculate GTE(T)
* using GTE(T) = [\sum_E GTE(E,T)*g(E)*e^{-E\Beta}] / [\sum_E g(E)*e^{-E\Beta}]
*
* Note both of the e terms has a hidden scaling constant (c) subtracted equal to
* max(log[g(e)] - E\Beta), such that the equations can be calculated in reasonable
* precision.
*
* In this version, we maintain 1024 histograms.
* We run one of more glauber cycles (with or without WL).
* For L=32, for each step of the glauber update, we accumulate in the e-th histogram
* (i.e. for steps (t, t+1), with energy values (e,e+1) which make up X,Y,W data, we store in hist e)
* For L > 32, we work out the hist as
* hist = { 0 if e = 0
*        { (e-1)/((L/32)^2) + 1 otherwise
*/
int sim_dos_gte_et(int argc, char* argv[])
{
	size_t L = 32, U = 1000, S = 1000, seed = 0;
	int runID = 0, reps = 1;
	double petThreshold = 0.0001, T = 1.0;
	std::string outputDir = "", dosFile = "";
	bool useSwendsenWang = false;

	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(outputDir, "directory")["--out-dir"]("Output directory")
			| Opt(L, "lattice")["-L"]["--lattice-size"]("Width and height of lattice")
			| Opt(S, "timesteps")["-S"]["--skip-steps"]("Number of skip steps")
			| Opt(U, "timesteps")["-U"]["--update-steps"]("Number of update steps")
			| Opt(seed, "seed")["--seed"]("Random seed (0 for unpredictable)")
			| Opt(runID, "ID")["--run-ID"]("ID used for labelling files")
			| Opt(dosFile, "filename")["--DoS-file"]("Input file containing density of states data")
			| Opt(T, "temperature")["-T"]["--temperature"]("Temperature value")
			| Opt(reps, "repetitions")["--reps"]("Ensemble repetitions")
			| Opt(useSwendsenWang)["--swendsen-wang"]("Use Swendsen-Wang updating for skip steps")
			| Opt(petThreshold, "threshold")["--Pet-threshold"]("Minimum threshold for P(E,T) values");
		auto result = cli.parse(Args(argc, argv));
		if (!result) {
			fmt::print("Error in command line: {}\n", result.errorMessage());
			return EXIT_FAILURE;
		}
		if (showHelp) {
			std::cerr << cli << std::endl;
			return EXIT_SUCCESS;
		}
	}
	
	fmt::print("Calculating per-glauber-step Tgl(E,T)\n");
	print_debug_info();

	size_t const N = L * L;
	std::vector<u8> lattice(N);
	std::vector<u8> lattice_buffer(N); // Used for accumulating G/TE hists

	//size_t const numE1 = 2 * L*L + 1;// 2 * 32 * 32 + 1;
	size_t const numE1 = 2 * 32 * 32 + 1;
	std::vector<Hist> gte_hists(numE1, Hist(numStates*numStates*numSiteEnergy, 0));
	Hist gte_counts(numE1, 0);
	std::vector<Hist> mi_hists(numE1, Hist(numStates*numStates, 0));
	Hist mi_counts(numE1, 0);
	std::vector<double> magnetisation(numE1, 0.);
	Hist mag_count(numE1, 0);
	std::vector<Hist> gte_binary_hists(numE1, Hist(numStates*numStates*16, 0));
	Hist gte_binary_counts(numE1, 0);
	std::vector<Hist> gte_hists_sweep(numE1, Hist(numStates*numStates*numSiteEnergy, 0));
	Hist gte_counts_sweep(numE1, 0);

	double invN = 1.0 / static_cast<double>(N);
	fmt::print("N={}, invN={}\n", N, invN);
	//size_t const numE = 2 * N + 1;
	
	std::vector<double> g;
	std::vector<int> e;
	size_t dosL;
	read_dos_file(dosFile.c_str(), g, e, dosL);
	if(dosL != L && dosL != 32) {
		PEEXIT("DoS File must use either use same lattice size or L=32\n");
	}
	fmt::print("dos file read\n");
	std::vector<double> pet(g.size());
	double maxVal = 0;
	for(int i = 0, len = static_cast<int>(g.size()); i < len; ++i) {
		pet[i] = g[i] - static_cast<double>(e[i]) / T;
		maxVal = std::max(maxVal, pet[i]);
	}
	for(auto &v : pet) {
		v = std::exp(v - maxVal);
	}
	pet.reserve(pet.size() + 4);
	//Fill in missing e values
	pet.insert(pet.end() - 2, 0);
	pet.insert(pet.end() - 1, 0);
	pet.insert(pet.end() - 1, 0);
	pet.insert(pet.end() - 1, 0);
	fmt::print("pet size: {}\n", pet.size());

	std::vector<bool> petMask(pet.size(), 0);
	for(int i = 0, len = static_cast<int>(petMask.size()); i < len; ++i) {
		petMask[i] = (pet[i] >= petThreshold);
	}

	auto usedSeed = seed ? seed : std::random_device()();
	std::mt19937_64 rand_engine(usedSeed);
	
	ConnectedSets sets(useSwendsenWang ? N : 0);

	std::uniform_real_distribution<> rng(0, 1);
	for(int rep = 0; rep < reps; ++rep) {

		initialise(rand_engine, lattice, L, (T < critical_temp('a')) ? 2: (rep % 2));

		//WangLandau wl(std::exp(1.0), std::exp(1e-08), 0.8, 1000, L, std::log(numStates), L*L * 2);
		int current_energy = calc_action(lattice, L);
		int prev_energy = 0;
		auto tpTable = tp_table_init(T, 'a');
		for (auto s = 0; s < S; ++s) {
			if (useSwendsenWang) {
				update_swendsen_wang(rand_engine, rng, lattice, L, N, T, current_energy, sets);
			}
			else {
				update_glauber(rand_engine, rng, lattice, L, N, tpTable, current_energy);
			}
		}
		for(auto u = 0; u < U; ++u) {
			lattice_buffer = lattice;
			prev_energy = current_energy;
			update_glauber_and_gte_hist(rand_engine, rng, lattice, L, N, tpTable, current_energy, petMask, gte_hists, gte_counts, mi_hists, mi_counts, magnetisation, mag_count, gte_binary_hists, gte_binary_counts, (dosL != L));
			GTE_reduced_histogram(lattice, lattice_buffer, L, gte_hists_sweep[scaleEnergyValue(prev_energy, L)]);
			if(u == U - 1) {
				int test_energy = calc_action(lattice, L);
				if(test_energy != current_energy) {
					fmt::print("Energy values have diverged: Running: {}, fresh: {}\n", current_energy, test_energy);
				}
			}
			auto scaledEnergy = scaleEnergyValue(current_energy, L);
			if(petMask[scaledEnergy]) {
				//magnetisation[scaledEnergy] += calc_magnetisation_proper(lattice, L);
				//++mag_count[scaledEnergy];
				for(auto y = 0u; y < L; ++y) {
					for(auto x = 0u; x < L; ++x) {
						int64_t s_i = lattice[y*L + x];
						int64_t s_r = lattice[y*L + wrap_plus(x, L - 1)];
						int64_t s_d = lattice[wrap_plus(y, L - 1)*L + x];
						int64_t s_l = lattice[y*L + wrap_minus(x, L - 1)];
						int64_t s_u = lattice[wrap_minus(y, L - 1)*L + x];

						++mi_hists[scaledEnergy][s_i + numStates * s_r];
						++mi_hists[scaledEnergy][s_i + numStates * s_l];
						++mi_hists[scaledEnergy][s_i + numStates * s_d];
						++mi_hists[scaledEnergy][s_i + numStates * s_u];
						mi_counts[scaledEnergy] += 4;
					}
				}
			}
		}
		progrep("Reps", rep, reps);
	}
	std::vector<double> gte_values(numE1, 0.);

	size_t totalCounts = 0;
	for(size_t i = 0, n = gte_hists.size(); i < n; ++i) {
		fmt::print("Calculating gte for hist {}, count: {}, petMask: {}\n", i, gte_counts[i], petMask[i] ? "true" : "false");
		gte_values[i] = calc_GTE_reduced_from_hist(gte_hists[i]);
		totalCounts += gte_counts[i];
	}
	fmt::print("Count: {} out of {}\n", totalCounts, L*L*U*reps);
	fmt::print("Done gte calculations\n");

	std::vector<double> gte_values_sweep(numE1, 0.);

	totalCounts = 0;
	for(size_t i = 0, n = gte_hists_sweep.size(); i < n; ++i) {
		fmt::print("Calculating gte for hist {}, count: {}, petMask: {}\n", i, gte_counts_sweep[i], petMask[i] ? "true" : "false");
		gte_values_sweep[i] = calc_GTE_reduced_from_hist(gte_hists_sweep[i]);
		totalCounts += gte_counts_sweep[i];
	}
	fmt::print("Count: {} out of {}\n", totalCounts, L*L*U*reps);
	fmt::print("Done gte_sweep calculations\n");

	std::vector<double> gte_binary_values(numE1, 0.);
	totalCounts = 0;
	for(size_t i = 0, n = gte_binary_hists.size(); i < n; ++i) {
		fmt::print("Calculating gte binary for hist {}, count: {}, petMask: {}\n", i, gte_binary_counts[i], petMask[i] ? "true" : "false");
		gte_binary_values[i] = calc_GTE_binary_from_hist(gte_binary_hists[i]);
		totalCounts += gte_binary_counts[i];
	}
	fmt::print("Count: {} out of {}\n", totalCounts, L*L*U*reps);
	fmt::print("Done gte calculations\n");

	totalCounts = 0;
	std::vector<double> mi_values(numE1, 0.);
	for(size_t i = 0, n = mi_hists.size(); i < n; ++i) {
		fmt::print("Calculating mi for hist {}, count: {}, petMask: {}\n", i, mi_counts[i], petMask[i] ? "true" : "false");
		mi_values[i] = calc_MI(mi_hists[i]);
		totalCounts += mi_counts[i];
	}
	fmt::print("Count: {} out of {}\n", totalCounts, L*L*U*reps);
	fmt::print("Done mi calculations\n");

	for(size_t i = 0; i < magnetisation.size(); ++i) {
		if(mag_count[i])
			magnetisation[i] /= mag_count[i];
	}

	std::string nameBuffer = fmt::format("{}/dosgteet_L{}_q{}_T{:.5}_{:04}.bin", outputDir, L, numStates, T, runID);
	FILE* binFile = fopen(nameBuffer.c_str(), "w");
	if(binFile == nullptr) PEEXIT("Failed to open binary file");

	magic_header(binFile, "potts stats");

	BIN_WRITE(CURRENT_GTE_FILE_VER, int, 1, binFile);
	BIN_WRITE(L, size_t, 1, binFile);
	BIN_WRITE(U, size_t, 1, binFile);
	BIN_WRITE(T, double, 1, binFile);
	BIN_WRITE(seed, size_t, 1, binFile);
	BIN_WRITE(numStates, u8, 1, binFile);
	BIN_WRITE(reps, int, 1, binFile);

	BIN_WRITE_PTR(&(gte_values[0]), double, numE1, binFile);
	BIN_WRITE_PTR(&(gte_counts[0]), size_t, numE1, binFile);
	BIN_WRITE_PTR(&(mi_values[0]), double, numE1, binFile);
	BIN_WRITE_PTR(&(mi_counts[0]), size_t, numE1, binFile);
	BIN_WRITE_PTR(&(magnetisation[0]), double, numE1, binFile);
	BIN_WRITE_PTR(&(mag_count[0]), size_t, numE1, binFile);
	BIN_WRITE_PTR(&(gte_binary_values[0]), double, numE1, binFile);
	BIN_WRITE_PTR(&(gte_binary_counts[0]), size_t, numE1, binFile);
	BIN_WRITE_PTR(&(gte_values_sweep[0]), double, numE1, binFile);
	BIN_WRITE_PTR(&(gte_counts_sweep[0]), size_t, numE1, binFile);

	fclose(binFile);

	return EXIT_SUCCESS;
}