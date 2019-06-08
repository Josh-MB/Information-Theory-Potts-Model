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
#include <fstream>

struct thread_context_sample
{
	std::vector<u8> lattice;
	std::vector<u8> lattice_buffer;
};

int sim_glauber(int argc, char* argv[])
{
	size_t L = 32, U = 1000, S = 1000, seed = 0, samples = 10;
	int runID = 0, threadLimit = 1, imode = 0, Tnum = 60, regime = 0;
	double Tmin = 0.001, Tmax = 1.7, orderedProportion = 1.0;
	char pottsVersion = 'a';
	std::string outputDir = "", TvalsFromFile = "";
	bool sortQStates = false, doRecordInterfaceLengths = false;

	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(outputDir, "directory")["--out-dir"]("Output directory")
			| Opt(L, "lattice")["-L"]["--lattice-size"]("Width and height of lattice")
			| Opt(U, "timesteps")["-U"]["--update-steps"]("Number of update steps")
			| Opt(S, "timesteps")["-S"]["--skip-steps"]("Number of skip steps")
			| Opt(seed, "seed")["--seed"]("Random seed (0 for unpredictable)")
			| Opt(runID, "ID")["--run-ID"]("ID used for labelling files")
			| Opt(threadLimit, "number")["--threads"]("Max number of threads to use")
			| Opt(samples, "samples")["--samples"]("Samples to take per temperature")
			| Opt(imode, "mode")["--init-mode"]("Initialisation mode")
			| Opt(Tnum, "count")["--T-count"]("How many temperature values to use")
			| Opt(Tmin, "temperature")["--T-min"]("Minimum temperature bound")
			| Opt(Tmax, "temperature")["--T-max"]("Maximum temperature bound")
			| Opt(pottsVersion, "type")["--potts-version"]("Potts version to use. Use b for Ising")
			| Opt(regime, "mode")["--heatbath-regime"]("Temperature regime. 0 - independent, 1 - warming, 2 - cooling")
			| Opt(orderedProportion, "threshold")["--low-T-ordered-prop"]("Proportion of lattice sites that are initialised to q=0 below T_c, if --init-mode=3")
			| Opt(sortQStates)["--sort-Q-states"]("Sorts states. Intended to overcome repetitions collapsing to different ground states")
			| Opt(TvalsFromFile, "filename")["--T-values"]("File to read temperature values from. Takes precedence over Tnum/Tmin/Tmax")
			| Opt(doRecordInterfaceLengths)["--interface-length"]("Calculate and record interface lengths");
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
	fmt::print("{:<12} = {:<24}\n", "--seed", seed);
	fmt::print("{:<12} = {:<24}\n", "--run-ID", runID);
	fmt::print("{:<12} = {:<24}\n", "--threads", threadLimit);
	fmt::print("{:<12} = {:<24}\n", "--samples", samples);
	fmt::print("{:<12} = {:<24}\n", "--init-mode", imode);
	fmt::print("{:<12} = {:<24}\n", "--T-count", Tnum);
	fmt::print("{:<12} = {:<24}\n", "--T-min", Tmin);
	fmt::print("{:<12} = {:<24}\n", "--T-max", Tmax);
	fmt::print("{:<12} = {:<24}\n", "--potts-version", pottsVersion);
	fmt::print("{:<12} = {:<24}\n", "--heatbath-regime", regime);
	fmt::print("{:<12} = {:<24}\n", "--low-T-ordered-prop", orderedProportion);
	fmt::print("{:<12} = {:<24}\n", "--sort-Q-states", sortQStates);
	fmt::print("{:<12} = {:<24}\n", "--T-values", TvalsFromFile);
	fmt::print("{:<12} = {:<24}\n", "--interface-length", doRecordInterfaceLengths);

	print_debug_info();

	size_t const N = L*L;

	if(pottsVersion != 'a' && pottsVersion != 'b') { ERRPT; fmt::print(stderr, "pottsVersion must be one of 'a' or 'b'"); exit(EXIT_FAILURE); }

	int num_threads = set_num_threads(threadLimit);

	// Prepare buffers for each thread
	std::vector<thread_context_sample> threadContexts(num_threads);
	for(auto &tc : threadContexts) {
		tc.lattice.resize(N);
		tc.lattice_buffer.resize(N);
		//prep_context_sample(tc);
	}

	// Calculate temperature values to use
	double T_crit = critical_temp(pottsVersion);
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
	// Warming
	if(regime == 1)
		std::sort(T_vals.begin(), T_vals.end());
	else if(regime == 2) // Cooling
		std::sort(T_vals.begin(), T_vals.end(), std::greater<double>());

	// Buffers for each temperature value
	std::vector<Hist> mi_hist(TnumActual, Hist(numStates*numStates, 0));
	std::vector<Hist> te_hist(TnumActual, Hist(numStates*numStates*numStates, 0));
	std::vector<Hist> gte_hist(TnumActual, Hist(static_cast<int64_t>(std::pow(numStates, 6)), 0));
	std::vector<Hist> gte_binary_hist(TnumActual, Hist(numStates* numStates* static_cast<int64_t>(std::pow(2, 4)), 0));
	std::vector<Hist> gte_reduced_hist(TnumActual, Hist(numStates*numStates*numSiteEnergy, 0));
	std::vector<Hist> energy_hist(TnumActual, Hist(2*N, 0));

	// Temporary buffers for each thread, gets accumulated into the above
	std::vector<Hist> working_mi_hist(num_threads, Hist(numStates* numStates, 0));
	std::vector<Hist> working_te_hist(num_threads, Hist(numStates* numStates* numStates, 0));
	std::vector<Hist> working_gte_hist(num_threads, Hist(static_cast<int64_t>(std::pow(numStates, 6)), 0));
	std::vector<Hist> working_energy_hist(num_threads, Hist(2 * N, 0));
	std::vector<Hist> working_gte_binary_hist(num_threads, Hist(numStates* numStates* static_cast<int64_t>(std::pow(2, 4)), 0));
	std::vector<Hist> working_gte_reduced_hist(num_threads, Hist(numStates*numStates*numSiteEnergy, 0));

	// Other buffers for magnetisation, autocorrelation time and interface length
	std::vector<double> avg_magnetisation(TnumActual);
	std::vector<std::vector<double>> working_magnetisation(num_threads);
	
	std::vector<double> avg_correlation_time(TnumActual);
	std::vector<double> avg_correlation_time_std_error(TnumActual);

	std::vector<double> avg_interfacial_lengths(TnumActual);
	std::vector<std::vector<int>> working_interfacial_lengths(num_threads);

	fmt::print("Initialised. Running simulations now\n");
	// Run T values independently, for a given number of samples
	// Parallel over T, not samples. This reduces how many threading
	// buffers are required i.e. things grouped according to temperature
	// like the histograms don't need thread buffers as only one thread
	// ever touches a single temperature value. 
	for(auto Ti = 0; Ti < TnumActual; ++Ti) {
		fmt::print("Running {}/{}, T={:.3}\n", Ti, TnumActual, T_vals[Ti]);
		
		// Initialise the transition probability look up table
		auto tpTable = tp_table_init(T_vals[Ti], pottsVersion);

		// Run samples
		#pragma omp parallel for
		for(auto sample = 0u; sample < samples; ++sample)
		{
			int const thread_num = omp_get_thread_num();
			auto& context = threadContexts[thread_num];

			std::mt19937_64 rand_engine(seed ? seed : std::random_device()());
			std::uniform_real_distribution<> rng(0, 1);

			// If below critical noise, start off magnetized
			int init = imode;
			if(imode == 3 || imode == 5 || imode == 6) {
				init = ((T_vals[Ti] < T_crit) ? 2 : 0);
				if((imode == 5) && (T_vals[Ti] >= T_crit) && (sample < samples / 2)) {
					init = 1;
					//init = 2;
				}
				if(imode == 6) {
					if((T_vals[Ti] >= T_crit - 0.00001 && T_vals[Ti] <= T_crit + 0.00001) && (sample < samples / 2)) {
						init = 2;
					}
				}
			}
			// Use ordered proportion value if below Tc
			if (imode == 7) {
				init = ((T_vals[Ti] < T_crit) ? 5 : 0);
			}
			// Only initialise if independent regime or if first temp value, or if cooling and below
			// critical temperature
			if(regime == 0 || Ti == 0 || (regime == 2 && T_vals[Ti] < T_crit)) {
				//printf("Initialising lattice (T=%f) with imode: %d\n", T_vals[Ti], init);
				initialise(rand_engine, context.lattice, L, init, orderedProportion);
			}

			// Run skip steps
			int current_energy = calc_action(context.lattice, L);
			for(size_t s = 0; s < S; ++s) {
				update_glauber(rand_engine, rng, context.lattice, L, N, tpTable, current_energy);
			}

			// Run update steps
			for(size_t u = 0; u < U; ++u) {
				// Save current state
				context.lattice_buffer = context.lattice;
				update_glauber(rand_engine, rng, context.lattice, L, N, tpTable, current_energy);
				// Accumulate data
				working_energy_hist[thread_num][current_energy]++;
				working_magnetisation[thread_num].emplace_back(calc_magnetisation(context.lattice, L));
				bool sortStates = (sortQStates == 1);// (T_vals[Ti] >= T_crit && current_energy > (1.25 * 2 * N));
				all_histograms(context.lattice, context.lattice_buffer, L,
							   working_mi_hist[thread_num], working_te_hist[thread_num], working_gte_hist[thread_num], working_gte_binary_hist[thread_num], working_gte_reduced_hist[thread_num], sortStates, rand_engine);
				if(doRecordInterfaceLengths) {
					working_interfacial_lengths[thread_num] = calcInterfaceAllLengths(context.lattice, static_cast<int>(L));
				}
			}
		}

		// Add all working histograms into the hist for this temperature
		// Also clear the working hists for the next temperature
		acc_thread_hists(working_mi_hist, mi_hist[Ti]);
		//acc_thread_hists(working_te_hist, te_hist[Ti]);
		acc_thread_hists(working_gte_hist, gte_hist[Ti]);
		//acc_thread_hists(working_gte_binary_hist, gte_binary_hist[Ti]);
		acc_thread_hists(working_gte_reduced_hist, gte_reduced_hist[Ti]);
		acc_thread_hists(working_energy_hist, energy_hist[Ti]);

		if(doRecordInterfaceLengths) {
			double length_total = 0.;
			size_t length_num = 0;
			for(auto& v : working_interfacial_lengths) {
				length_total += std::accumulate(v.begin(), v.end(), 0.);
				length_num += v.size();
				v.clear();
			}
			avg_interfacial_lengths[Ti] = length_total / length_num;
		}
		else {
			avg_interfacial_lengths[Ti] = -1.0;
		}

		double magnetisation_total = 0.;
		std::vector<double> tau;
		size_t num = 0;
		for(auto& v : working_magnetisation) {
			magnetisation_total += std::accumulate(v.begin(), v.end(), 0.);
			tau.emplace_back(calculate_relaxation_time(v));
			num += v.size();
			v.clear();
		}
		double tau_mean = std::accumulate(tau.begin(), tau.end(), 0.) / tau.size();
		double tau_var = std::accumulate(tau.begin(), tau.end(), 0.,
										 [tau_mean](double total, double x) {return total + (x - tau_mean)*(x - tau_mean); }) / tau.size();

		avg_magnetisation[Ti] = magnetisation_total / static_cast<double>(num);
		avg_correlation_time[Ti] = tau_mean;
		avg_correlation_time_std_error[Ti] = std::sqrt(tau_var) / tau.size();
	}

	std::string nameBuffer = fmt::format("{}/glauber_{:0.3}-{:0.3}_{:04}.bin", outputDir, Tmin, Tmax, runID);
	FILE* binFile = fopen(nameBuffer.c_str(), "w");
	if(binFile == nullptr) PEEXIT("Failed to open binary file");

	magic_header(binFile, "potts stats");

	BIN_WRITE(L, size_t, 1, binFile);
	BIN_WRITE(U, size_t, 1, binFile);
	BIN_WRITE(S, size_t, 1, binFile);
	BIN_WRITE(samples, size_t, 1, binFile);
	BIN_WRITE(imode, int, 1, binFile);
	BIN_WRITE(seed, size_t, 1, binFile);
	BIN_WRITE(numStates, u8, 1, binFile);
	BIN_WRITE(pottsVersion, char, 1, binFile);
	BIN_WRITE(regime, int, 1, binFile);
	BIN_WRITE(orderedProportion, double, 1, binFile);

	BIN_WRITE(TnumActual, int, 1, binFile);
	BIN_WRITE_PTR(&(T_vals[0]), double, static_cast<size_t>(TnumActual), binFile);

	//All histograms for all T values (i.e. TnumActual worth)
	for(auto& v: mi_hist)
		BIN_WRITE_PTR(&(v[0]), size_t, numStates*numStates, binFile);

	//for(auto& v : te_hist)
	//	BIN_WRITE_PTR(&(v[0]), size_t, numStates*numStates*numStates, binFile);

	for(auto& v : gte_hist)
		BIN_WRITE_PTR(&(v[0]), size_t, static_cast<size_t>(std::pow(numStates, 6)), binFile);

	//for(auto& v : gte_binary_hist)
	//	BIN_WRITE_PTR(&(v[0]), size_t, numStates*numStates*16, binFile);

	for(auto& v : gte_reduced_hist)
		BIN_WRITE_PTR(&(v[0]), size_t, numStates*numStates * numSiteEnergy, binFile);

	BIN_WRITE_PTR(&(avg_magnetisation[0]), double, static_cast<size_t>(TnumActual), binFile);
	BIN_WRITE_PTR(&(avg_interfacial_lengths[0]), double, static_cast<size_t>(TnumActual), binFile);

	for(auto& v : energy_hist)
		BIN_WRITE_PTR(&(v[0]), size_t, 2 * N, binFile);
	
	BIN_WRITE_PTR(&(avg_correlation_time[0]), double, static_cast<size_t>(TnumActual), binFile);
	BIN_WRITE_PTR(&(avg_correlation_time_std_error[0]), double, static_cast<size_t>(TnumActual), binFile);

	fclose(binFile);

	return EXIT_SUCCESS;
}
