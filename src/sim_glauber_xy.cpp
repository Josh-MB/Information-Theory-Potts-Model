#include "../include/model.hpp"
#include "../include/defs.hpp"
#include "../include/stats.hpp"
#include "../include/version.hpp"
#include "../include/utils.hpp"
#include "../include/dataSet.hpp"
#include "../include/psi.hpp"

#include <fmt/format.h>
#include <clara.hpp>
#include <vector>
#include <algorithm>
#include <random>
#include <omp.h>
#include <fstream>
#include <cstdio>

FILE* open_gnuplot_xy(char const*const gpcmd)
{
	FILE * gpipe = popen(gpcmd, "w");
	if(gpipe == nullptr) { ERRPT; fmt::print(stderr, "failed to open gpipe"); exit(EXIT_FAILURE); }
	if(setvbuf(gpipe, NULL, _IOLBF, 0) != 0) { ERRPT; fmt::print(stderr, "failed to line-buffer pipe to Gnuplot\n"); exit(EXIT_FAILURE); }
	return gpipe;
}

void init_frame_xy(FILE* gpipe, size_t const num_sub_plots, size_t const timestep, size_t const L, double const T, size_t const seed)
{
	std::string output_file;
	fmt::print(gpipe, "set term png size 1200,800\n");
	output_file = fmt::format("frames/xy_L{}_T{:.4f}_seed{}.png", L, T, seed);
	fmt::print(gpipe, "set out '{}'\n", output_file);

	if(num_sub_plots == 1) {
		fmt::print(gpipe, "set title 'MI for XY model: L={}, t={}\n", L, timestep);
	}
	else {
		fmt::print(gpipe, "set multiplot layout 1,{} title 'MI for XY model: L={}, t={}'\n", num_sub_plots, L, timestep);
	}
	fmt::print(gpipe, "unset key\n");
}

void finish_frame_xy(FILE* gpipe) 
{
	fmt::print(gpipe, "unset multiplot\n");
	fmt::print(gpipe, "unset out\n");
}

void draw_lattice(FILE* gpipe, size_t const L, std::vector<double> const& lattice)
{
	fmt::print(gpipe, "set title 'Lattice State'\n");
	fmt::print(gpipe, "set size square\n");
	fmt::print(gpipe, "set xr [-0.5:{}]\n", static_cast<double>(L) - 0.5);
	fmt::print(gpipe, "set yr [-0.5:{}]\n", static_cast<double>(L) - 0.5);
	fmt::print(gpipe, "unset xlabel\n");
	fmt::print(gpipe, "unset ylabel\n");
	fmt::print(gpipe, "plot '-' u ($2):($3):({0}*cos($1)):({0}*sin($1)) w vectors notitle\n", 0.1);
	for (size_t y = 0; y < L; ++y) {
		for (size_t x = 0; x < L; ++x) {
			fmt::print(gpipe, "{} {} {}\n", lattice[y*L+x], x, y);
		}
	}
	fprintf(gpipe, "e\n");
}

void draw_mi_xy(FILE* gpipe, size_t const L, const DataSet<2>& MI)
{
	VPoint<2> const* data = MI.getPoints();
	fmt::print(gpipe, "set title 'MI'\n");
	fmt::print(gpipe, "set size square\n");
	fmt::print(gpipe, "set xr [{}:{}]\n", -static_cast<double>(L), static_cast<double>(L));
	fmt::print(gpipe, "set yr [{}:{}]\n", -static_cast<double>(L), static_cast<double>(L));
	fmt::print(gpipe, "unset xlabel\n");
	fmt::print(gpipe, "unset ylabel\n");
	fmt::print(gpipe, "plot '-' u ($1):($2) w point\n");
	
	for (int i = 0; i < MI.getN(); ++i)
	{
		fmt::print(gpipe, "{} {}\n", data[i].values[0], data[i].values[1]);
	}
	fmt::print(gpipe, "e\n");
}

int sim_glauber_xy(int argc, char* argv[])
{
	size_t L = 32, U = 1000, S = 1000, seed = 0;
	int runID = 0, threadLimit = 1, Tnum = 60, estimatorNeighbours = 3;
	double Tmin = 0.001, Tmax = 1.7;
	bool calcMI = false, calcTE = false, calcGTE = false;
	std::string outputDir = "", TvalsFromFile = "";

	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(outputDir, "directory")["--out-dir"]("Output directory")
			| Opt(L, "lattice")["-L"]["--lattice-size"]("Width and height of lattice")
			| Opt(U, "timesteps")["-U"]["--update-steps"]("Number of update steps")
			| Opt(S, "timesteps")["-S"]["--skip-steps"]("Number of skip steps")
			| Opt(estimatorNeighbours, "neighbours")["-k"]["--estimator-neighbours"]("Number of neighbours to use in the continous entropy estimator")
			| Opt(seed, "seed")["--seed"]("Random seed (0 for unpredictable)")
			| Opt(runID, "ID")["--run-ID"]("ID used for labelling files")
			| Opt(threadLimit, "number")["--threads"]("Max number of threads to use")
			| Opt(Tnum, "count")["--T-count"]("How many temperature values to use")
			| Opt(Tmin, "temperature")["--T-min"]("Minimum temperature bound")
			| Opt(Tmax, "temperature")["--T-max"]("Maximum temperature bound")
			| Opt(calcMI)["--mi"]("Calculates MI metric")
			| Opt(calcTE)["--te"]("Calculates TE metric")
			| Opt(calcGTE)["--gte"]("Calculates GTE metric")
			| Opt(TvalsFromFile, "filename")["--T-values"]("File to read temperature values from. Takes precedence over Tnum/Tmin/Tmax");
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

	fmt::print("{:<12} = {:<24}\n", "--out-dir", outputDir);
	fmt::print("{:<12} = {:<24}\n", "--lattice-size", L);
	fmt::print("{:<12} = {:<24}\n", "--update-steps", U);
	fmt::print("{:<12} = {:<24}\n", "--skip-steps", S);
	fmt::print("{:<12} = {:<24}\n", "--estimator-neighbours", estimatorNeighbours);
	fmt::print("{:<12} = {:<24}\n", "--seed", seed);
	fmt::print("{:<12} = {:<24}\n", "--run-ID", runID);
	fmt::print("{:<12} = {:<24}\n", "--threads", threadLimit);
	fmt::print("{:<12} = {:<24}\n", "--T-count", Tnum);
	fmt::print("{:<12} = {:<24}\n", "--T-min", Tmin);
	fmt::print("{:<12} = {:<24}\n", "--T-max", Tmax);
	fmt::print("{:<12} = {:<24}\n", "--mi", calcMI);
	fmt::print("{:<12} = {:<24}\n", "--te", calcTE);
	fmt::print("{:<12} = {:<24}\n", "--gte", calcGTE);
	fmt::print("{:<12} = {:<24}\n", "--T-values", TvalsFromFile);

	print_debug_info();

	if (! (calcMI || calcTE || calcGTE)) {
		fmt::print("One of --mi, --te or --gte must be present\n");
		return EXIT_FAILURE;
	}
		
	int num_threads = set_num_threads(threadLimit);
	
	size_t const N = L*L;

	// Calculate temperature values to use
	// (2k_bT_c/J)ln(2k_bT_c/J) = 1
	// k_bT_c/J ~ 0.8816;
	constexpr double T_crit = 0.8816;
	std::vector<double> T_vals;
	if(!TvalsFromFile.empty()) {
		std::ifstream tfp(TvalsFromFile);
		double tval;
		while(tfp >> tval)
		{
			T_vals.push_back(tval);
		}
		tfp.close();
	}
	else
	{
		T_vals.resize(Tnum);
		if(Tmin > T_crit || Tmax < T_crit) { ERRPT; fmt::print(stderr, "bad temperature range (doesn't enclose critical temperature {})", T_crit); exit(EXIT_FAILURE); }
		double T_step = (Tmax - Tmin) / (Tnum - 1.0);
		double T = Tmin;
		int const i_crit = (int)round((T_crit - Tmin) / T_step);
		for(int64_t i = 0; i < Tnum; ++i, T += T_step) {
			T_vals[i] = T_crit + static_cast<double>(i - i_crit)*T_step;
		}
	}
	const int TnumActual = static_cast<int>(T_vals.size());

	std::vector<double> lattice(N, 0);
	std::vector<double> previousLattice(N, 0);
	std::vector<XY_site> latticeVectors;
	std::vector<double> siteEnergy(N, 0);
	ConnectedSets sets(N);
	
	// Since states might not flip, we'll end up with identical points
	// Because of this, we need to add low amplitude noise to break 
	// these degeneracies. See Kraskov, 2004, Phys Rev E
	// Unknown whether the same noise needs to be applied to points for
	// TE/GTE over multiple timesteps, so just in case, we'll save them
	constexpr int kDegeneracyNoiseAmplitude = -10;
	std::uniform_real_distribution<> degeneracy_rng(-pow(10.0, kDegeneracyNoiseAmplitude), pow(10.0, kDegeneracyNoiseAmplitude));
	std::vector<double> latticeDegeneracyBroken(N, 0);
	std::vector<double> previousLatticeDegeneracyBroken(N, 0);
	
	fmt::print("Initialised. Running simulations now\n");
		
	const int miN = calcMI ? N*U*2 : 0;
	const int teN = calcTE ? N*(U-1)*4 : 0;
	const int gteN = calcGTE ? N*(U-1) : 0;
	DataSet<2> dataForMI(miN, estimatorNeighbours);
	DataSet<3> dataForTE(teN, estimatorNeighbours);
	DataSet<6> dataForGTE(gteN, estimatorNeighbours);
	
	constexpr int Xdim = 0;
	constexpr int Ydim = 1;
	constexpr int WdimTE = 2;
	constexpr int WdimGTE = dataForGTE.dim - 1;
	
	std::vector<double> MI(TnumActual);
	std::vector<double> TE(TnumActual);
	std::vector<double> GTE(TnumActual);
	// Run T values independently, for a given number of samples
	// Parallel over T, not samples. This reduces how many threading
	// buffers are required i.e. things grouped according to temperature
	// like the histograms don't need thread buffers as only one thread
	// ever touches a single temperature value. 
	for(auto Ti = 0; Ti < TnumActual; ++Ti) {
		fmt::print("Running {}/{}, T={:.3}\n", Ti, TnumActual, T_vals[Ti]);
		
		std::mt19937_64 rand_engine(seed ? seed : std::random_device()());
		std::uniform_real_distribution<> rng(0, 1);

		std::uniform_real_distribution<> pi_rng(-PI, PI);
		for(auto &l : lattice) l = pi_rng(rand_engine);

		// Run skip steps
		for(size_t s = 0; s < S; ++s) {
			update_tomita(rand_engine, rng, lattice, L, N, T_vals[Ti], sets);
			progrep("Skip steps", s, S);
		}
		
		double current_energy = calc_energy_xy(lattice, latticeVectors, siteEnergy, L);

		// Run update steps
		int jMI = 0, jTE = 0, jGTE = 0;
		for(size_t u = 0; u < U; ++u) {
			// Save current state
			previousLattice = lattice;
			previousLatticeDegeneracyBroken = latticeDegeneracyBroken;
			update_glauber_xy(rand_engine, rng, lattice, latticeVectors, siteEnergy, T_vals[Ti], L, N, current_energy);

			progrep("Update steps", u, U);

			// Break degeneracies
			for (size_t i = 0; i < lattice.size(); ++i)
				latticeDegeneracyBroken[i] = lattice[i] + degeneracy_rng(rand_engine);
			const int L1 = L - 1;
			for (size_t y = 0; y < L; ++y) {
				for (size_t x = 0; x < L; ++x) {
					size_t idx = y*L+x;
					std::array<size_t, 4> neighbours{ wrap_minus(y, L1)*L + x, wrap_plus(y, L1)*L + x, y*L + wrap_minus(x, L1), y*L + wrap_plus(x, L1) };
					
					std::array<size_t, 2> mi_neighbours{ wrap_minus(y, L1)*L + x, y*L + wrap_minus(x, L1)};
					
					if (calcMI) {
						for (const size_t n : mi_neighbours) {
							dataForMI.setData(awrap(latticeDegeneracyBroken[idx]), jMI, Xdim);
							dataForMI.setData(awrap(latticeDegeneracyBroken[n]), jMI, Ydim);
							++jMI;
						}
					}
					
					if (calcTE) {
						for (const size_t n : neighbours) {
							if (u > 0) {
								dataForTE.setData(awrap(previousLatticeDegeneracyBroken[idx]), jTE, Xdim);
								dataForTE.setData(awrap(previousLatticeDegeneracyBroken[n]), jTE, Ydim);
								dataForTE.setData(awrap(latticeDegeneracyBroken[idx]), jTE, WdimTE);
								++jTE;
							}
						}
					}
					if (calcGTE && u > 0) {
						dataForGTE.setData(awrap(previousLatticeDegeneracyBroken[idx]), jGTE, Xdim);
						dataForGTE.setData(awrap(latticeDegeneracyBroken[idx]), jGTE, WdimGTE);
						for (int n = 0; n < neighbours.size(); ++n)
							dataForGTE.setData(awrap(previousLatticeDegeneracyBroken[neighbours[n]]), jGTE, Ydim + n);
						++jGTE;
					}
				}
			}
		}
		fmt::print("Added elements: {}, {}, {}\n", jMI, jTE, jGTE);
		
		/*const std::string gpcmd = "gnuplot";
		FILE* gpipe = open_gnuplot_xy(gpcmd.c_str());
		init_frame_xy(gpipe, 2, U+1, L, T_vals[Ti], seed);
		draw_lattice(gpipe, L, lattice);
		draw_mi_xy(gpipe, PI + 1.0, dataForMI);
		finish_frame_xy(gpipe);*/
		if (calcMI)	{
			fmt::print("Calculating MI\n");
			dataForMI.calcBallsizes();
			
			double psiNX = dataForMI.template countNeighboursOneDim<Xdim>();
			double psiNY = dataForMI.template countNeighboursOneDim<Ydim>();
			
			MI[Ti] = psi(dataForMI.getK()) + psi(dataForMI.getN()) - psiNX - psiNY;
			fmt::print("Done MI={} for T={}\n", MI[Ti], T_vals[Ti]);
		}
		if (calcTE)	{
			fmt::print("Calculating TE\n");
			dataForTE.calcBallsizes();
			
			double psiNXW = dataForTE.template countNeighboursTwoDim<Xdim, WdimTE>();
			double psiNXY = dataForTE.template countNeighboursTwoDim<Xdim, Ydim>();
			double psiNX = dataForTE.template countNeighboursOneDim<Xdim>();
			
			TE[Ti] = psi(dataForTE.getK()) - psiNXW - psiNXY + psiNX;
			fmt::print("Done TE={} for T={}\n", TE[Ti], T_vals[Ti]);
		}
		if (calcGTE) {
			fmt::print("Calculating GTE\n");
			dataForGTE.calcBallsizes();

			double psiNXW = dataForGTE.template countNeighboursTwoDim<Xdim, WdimGTE>();
			double psiNXY = dataForGTE.countNeighboursNDim();
			double psiNX = dataForGTE.template countNeighboursOneDim<Xdim>();

			GTE[Ti] = psi(dataForGTE.getK()) - psiNXW - psiNXY + psiNX;
			fmt::print("Done GTE={} for T={}\n", GTE[Ti], T_vals[Ti]);
		}
	}

	std::string nameBuffer = fmt::format("{}/glauber_xy_{:0.3}-{:0.3}_{:04}.bin", outputDir, Tmin, Tmax, runID);
	FILE* binFile = fopen(nameBuffer.c_str(), "wb");
	if(binFile == nullptr) PEEXIT("Failed to open binary file");

	magic_header(binFile, "xy stats");

	BIN_WRITE(L, size_t, 1, binFile);
	BIN_WRITE(U, size_t, 1, binFile);
	BIN_WRITE(S, size_t, 1, binFile);
	BIN_WRITE(seed, size_t, 1, binFile);
	BIN_WRITE(estimatorNeighbours,  int, 1, binFile);

	BIN_WRITE(TnumActual, int, 1, binFile);
	BIN_WRITE_PTR(&(T_vals[0]), double, static_cast<size_t>(TnumActual), binFile);

	BIN_WRITE_PTR(&(MI[0]), double, static_cast<size_t>(TnumActual), binFile);
	BIN_WRITE_PTR(&(TE[0]), double, static_cast<size_t>(TnumActual), binFile);
	BIN_WRITE_PTR(&(GTE[0]), double, static_cast<size_t>(TnumActual), binFile);
	fclose(binFile);

	return EXIT_SUCCESS;
}
