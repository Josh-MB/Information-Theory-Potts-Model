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
#include <bitset>

int sim_order(int argc, char* argv[])
{
	size_t L = 32, U = 1000, seed = 0;
	int runID = 0, reps = 1;
	std::string outputDir = "";
	double T = 1.0;

	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(outputDir, "directory")["--out-dir"]("Output directory")
			| Opt(L, "lattice")["-L"]["--lattice-size"]("Width and height of lattice")
			| Opt(U, "timesteps")["-U"]["--update-steps"]("Number of update steps")
			| Opt(seed, "seed")["--seed"]("Random seed (0 for unpredictable)")
			| Opt(runID, "ID")["--run-ID"]("ID used for labelling files")
			| Opt(reps, "repetitions")["--reps"]("Ensemble repetitions")
			| Opt(T, "temperature")["-T"]["--temp"]("Temperature to simulate at");
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

	constexpr int CURRENT_ORDER_FILE_VER = 1;

	fmt::print("Calculating order parameters\n");
	print_debug_info();

	size_t const N = L * L;
	std::vector<u8> lattice(N);

	auto usedSeed = seed ? seed : std::random_device()();
	std::mt19937_64 rand_engine(usedSeed);
	std::mt19937_64 second_rand_engine(usedSeed);
	std::uniform_real_distribution<> rng(0, 1);

	std::vector<double> magnetisation;
	auto numElements = reps * U;
	magnetisation.reserve(numElements);

	std::vector<double> one_mag;
	double tau_mean;
	double tau_std_error;
	one_mag.reserve(U);
	
	std::vector<double> tau;
	for(int rep = 0; rep < reps; ++rep) {
		initialise(rand_engine, lattice, L, (T < critical_temp('a')) ? 2 : (rep % 2));

		int current_energy = calc_action(lattice, L);
		auto tpTable = tp_table_init(T, 'a');

		for(auto u = 0u; u < U; ++u) {
			update_glauber(rand_engine, rng, lattice, L, N, tpTable, current_energy);
			auto mag = calc_magnetisation(lattice, L);
			magnetisation.emplace_back(mag);
			one_mag.emplace_back(mag);
		}
		tau.push_back(calculate_relaxation_time(one_mag));
		one_mag.clear();
	}
	tau_mean = std::accumulate(tau.begin(), tau.end(), 0.) / tau.size();
	double tau_var = std::accumulate(tau.begin(), tau.end(), 0.,
									 [tau_mean](double total, double x) {return total + (x - tau_mean)*(x - tau_mean); }) / tau.size();
	tau_std_error = std::sqrt(tau_var) / tau.size();

	double total_magnetisation = std::accumulate(magnetisation.begin(), magnetisation.end(), 0.);
	
	double avg_magnetisation = total_magnetisation / numElements;

	std::string nameBuffer = fmt::format("{}/gteorder_L{}_q{}_T{:.5}_{:04}.bin", outputDir, L, numStates, T, runID);
	FILE* binFile = fopen(nameBuffer.c_str(), "w");
	if(binFile == nullptr) PEEXIT("Failed to open binary file");

	magic_header(binFile, "potts order");

	BIN_WRITE(CURRENT_ORDER_FILE_VER, int, 1, binFile);
	BIN_WRITE(L, size_t, 1, binFile);
	BIN_WRITE(U, size_t, 1, binFile);
	BIN_WRITE(T, double, 1, binFile);
	BIN_WRITE(reps, int, 1, binFile);
	BIN_WRITE(seed, size_t, 1, binFile);
	BIN_WRITE(numStates, u8, 1, binFile);
	BIN_WRITE(avg_magnetisation, double, 1, binFile);
	BIN_WRITE(tau_mean, double, 1, binFile);
	BIN_WRITE(tau_std_error, double, 1, binFile);

	fclose(binFile);

	return EXIT_SUCCESS;
}