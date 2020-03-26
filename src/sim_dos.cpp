#include "../include/model.hpp"
#include "../include/defs.hpp"
#include "../include/stats.hpp"
#include "../include/wanglandau.hpp"
#include "../include/version.hpp"
#include "../include/utils.hpp"

#include <clara.hpp>
#include <fmt/format.h>

#include <vector>
#include <algorithm>
#include <random>
#include <omp.h>
#include <sstream>

/**
 * Used for calculating the density of states of the system
 */
int sim_dos(int argc, char* argv[])
{
	size_t L = 32, U = 1000, seed = 0;
	int runID = 0;
	std::string dosFile;
	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(L, "lattice")["-L"]["--lattice-size"]("Width and height of lattice")
			| Opt(U, "timesteps")["-U"]["--update-steps"]("Number of update steps")
			| Opt(seed, "seed")["--seed"]("Random seed (0 for unpredictable)")
			| Opt(runID, "ID")["--run-ID"]("ID used for labelling files")
			| Opt(dosFile, "filename")["--DoS-file"]("Input file containing density of states data for continuing a run");
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

	size_t const N = L*L;

	std::vector<u8> lattice(N);

	auto usedSeed = seed ? seed : std::random_device()();
	std::mt19937_64 rand_engine(usedSeed);
	
	std::uniform_real_distribution<> rng(0, 1);

	initialise(rand_engine, lattice, L, 0);

	WangLandau wl(std::exp(1.0), std::exp(1e-08), 0.8, U, L, std::log(numStates), L*L*2);

	// Load previous data
	if(!dosFile.empty()) {
		auto dosFileData = read_dos_file(dosFile.c_str(), &rand_engine);

		bool valid = true;
		CHECK_VAL(dosFileData.L, L, valid);
		CHECK_VAL(dosFileData.U, U, valid);
		if(seed != 0)
			CHECK_VAL(dosFileData.seed, seed, valid);
		if(!valid)
			PEEXIT("Input file is not valid");

		wl.setIteration(dosFileData.iteration);
		wl.setLogG(dosFileData.logG, dosFileData.factor);
	}

	int current_energy = calc_action(lattice, L);

	while(!wl.isDone()) {
		update_wang_landau(rand_engine, rng, lattice, wl, current_energy, L, N);
		if(wl.update()) {
			std::string nameBuffer = fmt::format("DoS_q{:d}_L{}_it{}_{:04d}.bin", numStates, L, wl.getIteration(), runID);
			FILE* binFile = fopen(nameBuffer.c_str(), "wb");
			if(binFile == nullptr) PEEXIT("Failed to open binary file");

			magic_header(binFile, "potts dos");

			BIN_WRITE(CURRENT_DOS_FILE_VER, int, 1, binFile);
			BIN_WRITE(L, size_t, 1, binFile);
			BIN_WRITE(U, size_t, 1, binFile);
			BIN_WRITE(usedSeed, size_t, 1, binFile);
			std::stringbuf sbuf;
			std::ostream os(&sbuf);
			os << rand_engine;

			auto len = strlen(sbuf.str().c_str());
			BIN_WRITE(len, size_t, 1, binFile);
			BIN_WRITE_PTR(sbuf.str().c_str(), char, len, binFile);
			BIN_WRITE(numStates, u8, 1, binFile);
			auto factor = wl.getFactor();
			BIN_WRITE(factor, double, 1, binFile);
			auto it = wl.getIteration();
			BIN_WRITE(it, size_t, 1, binFile);

			auto logG = wl.getLogG();
			auto finalHist = wl.getFinalHist();
			size_t numZeroes = 0u;
			for(auto h : finalHist) {
				if(h == 0) ++numZeroes;
			}

			auto fullSz = logG.size();
			size_t sz = fullSz - numZeroes;

			BIN_WRITE(sz, size_t, 1, binFile);
			BIN_WRITE_PTR(&(logG[0]), double, logG.size(), binFile);
			BIN_WRITE_PTR(&(finalHist[0]), size_t, finalHist.size(), binFile);
			
			for(size_t i = 0u; i < fullSz; ++i) {
				size_t j = i * 2;
				BIN_WRITE(j, size_t, 1, binFile);
			}

			BIN_WRITE_PTR(&(lattice[0]), u8, lattice.size(), binFile);
			//Need to record energy levels too
			//Only record g where h>0
			fclose(binFile);
		}
	}

	return EXIT_SUCCESS;
}