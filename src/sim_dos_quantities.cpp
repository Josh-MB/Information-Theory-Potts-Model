#include "../include/model.hpp"
#include "../include/defs.hpp"
#include "../include/stats.hpp"
#include "../include/wanglandau.hpp"
#include "../include/version.hpp"
#include "../include/utils.hpp"

#include <fmt/format.h>
#include <clara.hpp>
#include <vector>
#include <algorithm>
#include <random>
#include <sstream>
#include <omp.h>

/**
* Used for calculating the magnetisation and MI of the system
* using a WL update scheme. This is possible due to these quantities
* being instantaneous, and thus we can just measure at every E rather than
* needing the "true" sequence (i.e., that given by a glauber update)
*/
int sim_dos_quantities(int argc, char* argv[])
{
	size_t L = 32, U = 1000, seed = 0, histCount = 5;
	int runID = 0, reps = 1, threadLimit = 1, imode = 0;
	std::string outputDir = "", dosFile = "";
	bool useOldG = false;
	double sortThreshold = 1.0;

	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(outputDir, "directory")["--out-dir"]("Output directory")
			| Opt(L, "lattice")["-L"]["--lattice-size"]("Width and height of lattice")
			| Opt(U, "timesteps")["-U"]["--update-steps"]("Number of update steps")
			| Opt(seed, "seed")["--seed"]("Random seed (0 for unpredictable)")
			| Opt(runID, "ID")["--run-ID"]("ID used for labelling files")
			| Opt(dosFile, "filename")["--DoS-file"]("Initialise logG data with previous run")
			| Opt(threadLimit, "number")["--threads"]("Max number of threads to use")
			| Opt(reps, "repetitions")["--reps"]("Ensemble repetitions")
			| Opt(imode, "mode")["--init-mode"]("Initialisation mode")
			| Opt(useOldG)["--use-prev-logG"]("Use g weights in DoS file")
			| Opt(sortThreshold, "threshold")["--sort-threshold"]("Sorts states when magnetisation is above threshold. Intended to overcome low temperature repetitions collapsing to different ground states")
			| Opt(histCount, "count")["--E-visits"]("Update until every energy value has been visited this many times");
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

	int constexpr CURRENT_FILE_VER = 1;

	int num_threads = set_num_threads(threadLimit);

	print_debug_info();

	size_t const N = L * L;
	size_t const E = N * 2;
	size_t const E1 = E + 1;

	std::vector<u8> lattice(N);

	auto usedSeed = seed ? seed : std::random_device()();
	std::mt19937_64 rand_engine(usedSeed);

	std::uniform_real_distribution<> rng(0, 1);

	std::vector<u8> backupLattice(N);
	
	double factor1 = 1.0000000074505806;
	std::vector<double> logG(E1, 0);

	if(!dosFile.empty()) {
		auto dosFileData = read_dos_file(dosFile.c_str(), &rand_engine);
		
		bool valid = true;
		CHECK_VAL(dosFileData.L, L, valid);
		if (!valid)
			PEEXIT("Input file is not valid");

		backupLattice = dosFileData.lattice;
		factor1 = dosFileData.factor;
		logG = dosFileData.logG;
	}
	else {
		PEEXIT("No DoS file provided");
	}

	std::vector<Hist> mi_hists(E1, Hist(numStates*numStates, 0));
	std::vector<double> mag_values(E1, 0);
	Hist countsTotal(E1, 0);
	std::vector<double> inst_mi_vals(E1, 0);
	Hist interface_running_totals(E1, 0);
	Hist interface_running_num(E1, 0);

	std::vector<std::vector<u8>> latticeQueue(num_threads, std::vector<u8>(N, 0));
	std::vector<int> energyQueue(num_threads);

	for(int r = 0; r < reps; ++r) {
		Hist counts(E1, 0);
		// Ising model can't have odd energy, or energy = 2, or
		// energy = E-2, so ignore those states
		if (numStates == 2) {
			for(int e = 0; e < E1; ++e) {
				if(((e % 2) == 1) || (e == 2) || (e == (E - 2)))
				   counts[e] = histCount * 2;
			}
		}
		// Otherwise, q-state Potts (for q>2) can't have the following
		// energy levels
		else {
			counts[E1 - 2] = histCount * 2;
			counts[E1 - 3] = histCount * 2;
			counts[E1 - 4] = histCount * 2;
			counts[E1 - 6] = histCount * 2;
		}
		//counts[0] = histCount * 2;
		auto minCount = std::min_element(counts.begin(), counts.end());


		if(imode == 5) {
			lattice = backupLattice;
		}
		else {
			auto initVal = imode;
			if(imode == 3)
				initVal = (r % 2);
			initialise(rand_engine, lattice, L, initVal);
		}
		int current_energy = calc_action(lattice, L);

		WangLandau wl(std::exp(1.0), std::exp(1e-08), 0.8, U, L, std::log(numStates), E);
		if(useOldG) {
			wl.setLogG(logG, factor1);
		}

		size_t iteration = 0;
		while(*minCount < histCount) {
			bool updateMinCount = update_wang_landau_and_quantities(rand_engine, rng, lattice, wl, current_energy, L, N,
																	mi_hists, mag_values, counts, static_cast<int>(histCount * 2), inst_mi_vals, static_cast<int>(std::distance(counts.begin(), minCount)), sortThreshold,
																	interface_running_totals, interface_running_num, latticeQueue, num_threads, energyQueue);
			
			// If the position pointed to by the minCount iterator was updated
			// we need to find the new minimum
			if(updateMinCount) { //std::distance(counts.begin(), minCount) == current_energy) {
				minCount = std::min_element(counts.begin(), counts.end());
			}

			++iteration;
			if(iteration % 10000 == 0)
			{
				auto countSum = std::accumulate(counts.begin(), counts.end(), size_t{});
				if(countSum == E1 * histCount * 2) break;
				fmt::print("Iteration (new) {}, minCount: {} (element: {}), energy: {}, countSum: {}\n", iteration, *minCount, static_cast<int>(std::distance(counts.begin(), minCount)), current_energy, countSum);
			}
		}

		std::transform(counts.begin(), counts.end(), countsTotal.begin(), countsTotal.begin(), std::plus<size_t>());
	}
	std::string miNameBuffer = fmt::format("{}/mi_hists_L{}_q{}_{:04}.bin", outputDir, L, numStates, runID);
	FILE* miFile = fopen(miNameBuffer.c_str(), "w");
	if(miFile == nullptr) PEEXIT("Failed to open mi file");
	std::vector<double> mi_values(E1, 0.);
	for(size_t i = 0, n = mi_hists.size(); i < n; ++i) {
		mi_values[i] = calc_MI(mi_hists[i]);
		BIN_WRITE_PTR(&(mi_hists[i][0]), size_t, numStates*numStates, miFile);
	}
	fmt::print("Done mi calculations\n");
	fclose(miFile);

	for(size_t i = 0; i < mag_values.size(); ++i) {
		if(countsTotal[i] == 0)  PEEXIT("Count %d is 0. This shouldn't be possible", i);
		mag_values[i] /= countsTotal[i];
	}

	std::vector<double> avg_interface_lengths(E1, 0.);
	for(int e = 0; e < interface_running_totals.size(); ++e) {
		avg_interface_lengths[e] = static_cast<double>(interface_running_totals[e]) / interface_running_num[e];
	}

	/*for(size_t i = 0; i < inst_mi_vals.size(); ++i) {
		if(counts[i] == 0)  PEEXIT("Count %d is 0. This shouldn't be possible", i);
		inst_mi_vals[i] /= counts[i];
	}*/

	std::string nameBuffer = fmt::format("{}/dos_quantities_L{}_q{}_{:04}.bin", outputDir, L, numStates, runID);
	FILE* binFile = fopen(nameBuffer.c_str(), "w");
	if(binFile == nullptr) PEEXIT("Failed to open binary file");

	magic_header(binFile, "potts dos_quantities");

	BIN_WRITE(CURRENT_FILE_VER, int, 1, binFile);
	BIN_WRITE(L, size_t, 1, binFile);
	BIN_WRITE(U, size_t, 1, binFile);
	BIN_WRITE(seed, size_t, 1, binFile);
	BIN_WRITE(numStates, u8, 1, binFile);
	BIN_WRITE(histCount, int, 1, binFile);

	BIN_WRITE_PTR(&(mi_values[0]), double, E1, binFile);
	BIN_WRITE_PTR(&(mag_values[0]), double, E1, binFile);
	BIN_WRITE_PTR(&(countsTotal[0]), size_t, E1, binFile);
	BIN_WRITE_PTR(&(inst_mi_vals[0]), double, E1, binFile);
	BIN_WRITE_PTR(&(avg_interface_lengths[0]), double, E1, binFile);

	fclose(binFile);

	return EXIT_SUCCESS;
}