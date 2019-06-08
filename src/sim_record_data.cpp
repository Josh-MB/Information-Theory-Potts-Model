#include "../include/model.hpp"
#include "../include/defs.hpp"
#include "../include/stats.hpp"
#include "../include/version.hpp"
#include "../include/utils.hpp"

#include <fmt/format.h>
#include <clara.hpp>
#include <vector>
#include <algorithm>
#include <random>
#include <omp.h>

int sim_record_data(int argc, char* argv[])
{
	size_t L = 32, U = 10, S = 1000, seed = 0;
	int runID = 0, reps = 1, imode = 2;
	std::string outputDir = "";
	double T = 1.0;

	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(outputDir, "directory")["--out-dir"]("Output directory")
			| Opt(L, "lattice")["-L"]["--lattice-size"]("Width and height of lattice")
			| Opt(U, "timesteps")["-U"]["--update-steps"]("Number of update steps")
			| Opt(S, "timesteps")["-S"]["--skip-steps"]("Number of skip steps")
			| Opt(imode, "mode")["--init-mode"]("Initialisation mode")
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

	print_debug_info();

	size_t const N = L * L;

	std::vector<u8> lattice(N, 0);
	std::vector<u8> lattice_buffer(N, 0);

	auto tpTable = tp_table_init(T, 'a');

	std::mt19937_64 rand_engine(seed ? seed : std::random_device()());
	std::uniform_real_distribution<> rng(0, 1);

	std::string nameBuffer = fmt::format("{}/record_data_q{}_L{}_T{:0.3}_{:04}.bin", outputDir, numStates, L, T, runID);
	FILE* binFile = fopen(nameBuffer.c_str(), "w");
	if(binFile == nullptr) PEEXIT("Failed to open binary file");

	magic_header(binFile, "potts record");
	BIN_WRITE(L, size_t, 1, binFile);
	BIN_WRITE(U, size_t, 1, binFile);
	BIN_WRITE(S, size_t, 1, binFile);
	BIN_WRITE(imode, int, 1, binFile);
	BIN_WRITE(seed, size_t, 1, binFile);
	BIN_WRITE(numStates, u8, 1, binFile);
	BIN_WRITE(T, double, 1, binFile);
	BIN_WRITE(reps, int, 1, binFile);

	for(int r = 0; r < reps; ++r) {
		int init = (reps % 2) * 2;
		initialise(rand_engine, lattice, L, init);

		// Run skip steps
		int ignore_current_energy = 0;
		for(size_t s = 0; s < S; ++s) {
			update_glauber(rand_engine, rng, lattice, L, N, tpTable, ignore_current_energy);
		}

		// Run update steps
		for(size_t u = 0; u < U; ++u) {
			// Save current state
			fmt::print("u={}\n", u);
			lattice_buffer = lattice;
			BIN_WRITE_PTR(&(lattice[0]), u8, N, binFile);
			update_glauber_and_record(rand_engine, rng, lattice, L, N, tpTable, ignore_current_energy, binFile);

			size_t numSame = 0;
			for(int i = 0; i < lattice.size(); ++i) {
				numSame += (lattice[i] == lattice_buffer[i]);
			}
			fmt::print("Num same (per sweep): {}\n", numSame);
		}
	}
	fclose(binFile);
	return EXIT_SUCCESS;
}