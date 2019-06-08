#include "../include/model.hpp"
#include "../include/defs.hpp"
#include "../include/stats.hpp"
#include "../include/utils.hpp"

#include <fmt/format.h>
#include <clara.hpp>
#include <vector>
#include <algorithm>
#include <random>
#include <omp.h>

int sim_kl_filtered_gte(int argc, char* argv[])
{
	size_t L = 32, U = 1000, S = 1000, seed = 0, samples = 10, KLLag = 10, KLWindows = 100;
	int runID = 0, threadLimit = 1, imode = 0, Tnum = 60, regime = 0;
	double Tmin = 0.001, Tmax = 1.7, KLThreshold = 0.5;
	char pottsVersion = 'a';
	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
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
			| Opt(KLThreshold, "threshold")["--KL-threshold"]("Threshold for GTE Kullback-Leibler Divergence filtering")
			| Opt(KLWindows, "count")["--KL-windows"]("Number of KL windows to test and form filtered histograms from")
			| Opt(KLLag, "count")["--KL-lag"]("Window lag to compare against for KL divergence (i.e D_KL(hist(t-lag)||hist(t))");
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
	double const inv_L_sq = 1.0 / static_cast<double>(N);

	if(pottsVersion != 'a' && pottsVersion != 'b') { ERRPT; fmt::print(stderr, "pottsVersion must be one of 'a' or 'b'"); exit(EXIT_FAILURE); }

	std::vector<u8> lattice(N);
	std::vector<u8> lattice_buffer(N);

	// Calculate temperature values to use
	double T_crit = critical_temp(pottsVersion);
	if(Tmin > T_crit || Tmax < T_crit) { ERRPT; fmt::print(stderr, "bad temperature range (doesn't enclose critical temperature {})", T_crit); exit(EXIT_FAILURE); }
	double T_step = (Tmax - Tmin) / (Tnum - 1.0);
	std::vector<double> T_vals(Tnum);
	double T = Tmin;
	int const i_crit = (int)round((T_crit - Tmin) / T_step);
	for(int64_t i = 0; i < Tnum; ++i, T += T_step) {
		T_vals[i] = T_crit + static_cast<double>(i - i_crit)*T_step;
	}

	// Warming
	if(regime == 1)
		std::sort(T_vals.begin(), T_vals.end());
	else if(regime == 2) // Cooling
		std::sort(T_vals.begin(), T_vals.end(), std::greater<double>());

	std::vector<Hist> mi_hist(Tnum, Hist(numStates * numStates, 0));
	std::vector<Hist> te_hist(Tnum, Hist(numStates * numStates * numStates, 0));
	std::vector<Hist> gte_hist(Tnum, Hist(static_cast<int64_t>(std::pow(numStates, 6)), 0));
	std::vector<Hist> gte_binary_hist(Tnum, Hist(numStates * numStates * 16, 0));
	std::vector<Hist> gte_reduced_hist(Tnum, Hist(numStates*numStates*numNeighbourhoods, 0));

	Hist working_mi_hist(numStates*numStates, 0);
	Hist working_te_hist(numStates*numStates*numStates, 0);
	Hist working_gte_hist(static_cast<int64_t>(std::pow(numStates, 6)), 0);
	Hist working_gte_binary_hist(numStates*numStates * 16, 0);
	Hist working_gte_reduced_hist(numStates*numStates * numNeighbourhoods, 0);

	std::vector<Hist> gte_hist_KL_lo(Tnum, Hist(static_cast<int64_t>(std::pow(numStates, 6)), 0));
	std::vector<Hist> gte_hist_KL_hi(Tnum, Hist(static_cast<int64_t>(std::pow(numStates, 6)), 0));
	std::vector<Hist> gte_hist_buffer(KLLag, Hist(static_cast<int64_t>(std::pow(numStates, 6)), 0));

	std::vector<double> working_energy;
	std::vector<int> working_action(U);

	std::vector<double> avg_energy(Tnum);
	std::vector<std::vector<int>> action(Tnum, std::vector<int>(samples * U, 0));

	for(auto Ti = 0; Ti < Tnum; ++Ti) {
		fmt::print("Running {}/{}, T={:.3}\n", Ti, Tnum, T_vals[Ti]);

		// Initialise the transition probability look up table
		auto tpTable = tp_table_init(T_vals[Ti], pottsVersion);
		for(auto sample = 0u; sample < samples; ++sample)
		{
			std::mt19937_64 rand_engine(seed ? seed : std::random_device()());
			std::uniform_real_distribution<> rng(0, 1);

			int init = imode;
			if(imode == 3) {
				// If below critical noise, start off magnetized
				init = ((T_vals[Ti] < T_crit) ? 2 : 0);
			}
			// Only initialise if independent regime or if first temp value, or if cooling and below
			// critical temperature
			if(regime == 0 || Ti == 0 || (regime == 2 && T_vals[Ti] < T_crit)) {
				//printf("Initialising lattice (T=%f) with imode: %d\n", T_vals[Ti], init);
				initialise(rand_engine, lattice, L, init);
			}

			// Run skip steps
			int ignore_current_energy = 0;
			for(size_t s = 0; s < S; ++s) {
				update_glauber(rand_engine, rng, lattice, L, N, tpTable, ignore_current_energy);
			}

			for(size_t kl = 0; kl < KLWindows; ++kl) {
				printf("kl=%zu\n", kl);
				// Run update steps
				for(size_t u = 0; u < U; ++u) {
					// Save current state
					lattice_buffer = lattice;
					update_glauber(rand_engine, rng, lattice, L, N, tpTable, ignore_current_energy);
					//lattice_history[u] = lattice;
					working_energy.emplace_back(calculate_total_energy(lattice, L, pottsVersion) * inv_L_sq);
					working_action[u] = calc_action(lattice, L);

					all_histograms(lattice, lattice_buffer, L,
								   working_mi_hist, working_te_hist, working_gte_hist, working_gte_binary_hist, working_gte_reduced_hist, false, rand_engine);
				}
				if(kl >= KLLag) {
					double kl_val = kl_divergence(gte_hist_buffer[0], working_gte_hist);
					Hist* gteFiltered = ((kl_val < KLThreshold) ? &(gte_hist_KL_lo[Ti]) : &(gte_hist_KL_hi[Ti]));
					acc_hist(working_gte_hist, *gteFiltered);
					//Update histogram buffer
					for(size_t i = 0u; i < KLLag - 1; ++i) {
						gte_hist_buffer[i] = gte_hist_buffer[i + 1];
					}
					gte_hist_buffer[KLLag - 1] = working_gte_hist;
				}
				else {
					gte_hist_buffer[kl] = working_gte_hist;
				}

				if(kl == KLWindows - 1) {
					acc_hist(working_mi_hist, mi_hist[Ti]);
					acc_hist(working_te_hist, te_hist[Ti]);
					acc_hist(working_gte_hist, gte_hist[Ti]);
					acc_hist(working_gte_binary_hist, gte_binary_hist[Ti]);

					avg_energy[Ti] = std::accumulate(working_energy.begin(), working_energy.end(), 0.);
					working_energy.clear();

					action[Ti].insert(action[Ti].end(), working_action.begin(), working_action.end());
				}
				clear_and_resize_vec(working_mi_hist);
				clear_and_resize_vec(working_te_hist);
				clear_and_resize_vec(working_gte_hist);
				clear_and_resize_vec(working_gte_binary_hist);
				clear_and_resize_vec(working_action);
			}
		}
	}

	std::transform(avg_energy.begin(), avg_energy.end(), avg_energy.begin(), [U, samples](double e) {return e / (U*samples); });

	std::string nameBuffer = fmt::format("glauberkl_{:0.3}-{:0.3}_{:04}.bin", Tmin, Tmax, runID);
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

	BIN_WRITE(Tnum, int, 1, binFile);
	BIN_WRITE_PTR(&(T_vals[0]), double, static_cast<size_t>(Tnum), binFile);

	BIN_WRITE_PTR(&(avg_energy[0]), double, static_cast<size_t>(Tnum), binFile);

	//All histograms for all T values (i.e. Tnum worth)
	for(auto& v : mi_hist)
		BIN_WRITE_PTR(&(v[0]), size_t, numStates*numStates, binFile);

	for(auto& v : te_hist)
		BIN_WRITE_PTR(&(v[0]), size_t, numStates*numStates*numStates, binFile);

	for(auto& v : gte_hist)
		BIN_WRITE_PTR(&(v[0]), size_t, static_cast<size_t>(std::pow(numStates, 6)), binFile);

	for(auto& v : action)
		BIN_WRITE_PTR(&(v[0]), int, U*samples, binFile);

	for(auto& v : gte_binary_hist)
		BIN_WRITE_PTR(&(v[0]), size_t, numStates*numStates * 16, binFile);

	for(auto& v : gte_hist_KL_lo)
		BIN_WRITE_PTR(&(v[0]), size_t, static_cast<size_t>(std::pow(numStates, 6)), binFile);

	for(auto& v : gte_hist_KL_hi)
		BIN_WRITE_PTR(&(v[0]), size_t, static_cast<size_t>(std::pow(numStates, 6)), binFile);

	fclose(binFile);

	return EXIT_SUCCESS;
}