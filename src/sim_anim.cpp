#include "../include/model.hpp"
#include "../include/defs.hpp"
#include "../include/stats.hpp"
#include "../include/utils.hpp"

#include <clara.hpp>
#include <fmt/format.h>
#include <vector>
#include <algorithm>
#include <random>

FILE* open_gnuplot(char const*const gpcmd)
{
	FILE * gpipe = popen(gpcmd, "w");
	if(gpipe == nullptr) { ERRPT; fmt::print(stderr, "failed to open gpipe"); exit(EXIT_FAILURE); }
	if(setvbuf(gpipe, NULL, _IOLBF, 0) != 0) { ERRPT; fmt::print(stderr, "failed to line-buffer pipe to Gnuplot\n"); exit(EXIT_FAILURE); }
	return gpipe;
}

void init_frame(FILE* gpipe, size_t const num_sub_plots, size_t const timestep, size_t const L, size_t const seed, int const filetype)
{
	std::string output_file;
	if(filetype == 1) {
		fmt::print(gpipe, "set terminal postscript eps enhanced size 6in,6in\n");
		output_file = fmt::format("frames/potts_{}_seed{}.ps", timestep, seed);
	}
	else if(filetype == 0) {
		fmt::print(gpipe, "set term png size 1200,800\n");
		output_file = fmt::format("frames/potts_{}_seed{}.png", timestep, seed);
	}
	fmt::print(gpipe, "set out '{}'\n", output_file);

	if(num_sub_plots == 1) {
		fmt::print(gpipe, "set title '5-state Potts Model: L={}, t={}'\n", L, timestep);
	}
	else {
		fmt::print(gpipe, "set multiplot layout 1,{} title '5-state Potts Model: L={}, t={}'\n", num_sub_plots, L, timestep);
	}
	fmt::print(gpipe, "unset key\n");
}
void finish_frame(FILE* gpipe) 
{
	fmt::print(gpipe, "unset multiplot\n");
	fmt::print(gpipe, "unset out\n");
}

void draw_lattice(FILE* gpipe, size_t const L, std::vector<u8> const& lattice)
{
	fmt::print(gpipe, "set title 'Lattice State'\n");
	fmt::print(gpipe, "set size square\n");
	fmt::print(gpipe, "set xr [-0.5:{}]\n", static_cast<double>(L) - 0.5);
	fmt::print(gpipe, "set yr [-0.5:{}]\n", static_cast<double>(L) - 0.5);
	fmt::print(gpipe, "unset xlabel\n");
	fmt::print(gpipe, "unset ylabel\n");
	fmt::print(gpipe, "set palette model RGB maxcolors 5\n");
	fmt::print(gpipe, "set palette defined ( 0 '#E31A1C', 1 '#FF7F00', 2 '#33A02C', 3 '#1F78B4', 4 '#6A3D9A', 5 '#FB9A99', 6 '#FDBF6F', 7 '#B2DF8A', 8 '#A6CEE3', 9 '#CAB2D6')\n");
	fmt::print(gpipe, "unset colorbox\n");
	fmt::print(gpipe, "plot '-' u 1:2:3:4:5 with rgbimage notitle\n");

	int const reds[] =		{ 207, 255, 51,  31,  106, 251, 253, 178, 166, 202};
	int const greens[] =	{ 26,  127, 160, 120, 61,  154, 191, 223, 206, 178};
	int const blues[] =		{ 28,  0,   44,  180, 154, 153, 111, 138, 227, 214};

	for(size_t y = 0; y < L; ++y) {
		for(size_t x = 0; x < L; ++x) {
			size_t state = static_cast<size_t>(lattice[y*L + x]);
			fmt::print(gpipe, "{} {} {} {} {}\n", x, y, reds[state], greens[state], blues[state]);
		}
		fmt::print(gpipe, "\n");
	}
	fmt::print(gpipe, "\ne\n");
}
void draw_energy(FILE* gpipe, size_t const max_time_step, std::vector<double> const& energy, std::vector<double> const& mi)
{
	fmt::print(gpipe, "set title 'Total Energy: {}, MI: {}'\n", energy.back(), mi.back());
	fmt::print(gpipe, "set size square\n");
	fmt::print(gpipe, "set xr [0:{}]\n", max_time_step);
	fmt::print(gpipe, "set yr [-4.5:1]\n");
	fmt::print(gpipe, "set border\n");
	fmt::print(gpipe, "set tics\n");
	fmt::print(gpipe, "set key\n");
	fmt::print(gpipe, "set xlabel 'Time step'\n");
	fmt::print(gpipe, "set ylabel 'Total Energy'\n");
	fmt::print(gpipe, "plot '-' u ($1):($2) w lines title 'Energy' lc rgb 'red', '-' u ($1):($2) w lines title 'Instantaneous MI' lc rgb 'blue'\n");

	auto sz = energy.size();
	for(auto x = 0u; x < sz; ++x) {
		fmt::print(gpipe, "{} {}\n", x, energy[x]);
	}
	fmt::print(gpipe, "e\n");

	sz = mi.size();
	for(auto x = 0u; x < sz; ++x) {
		fmt::print(gpipe, "{} {}\n", x, mi[x]);
	}
	fmt::print(gpipe, "e\n");
}

int sim_anim(int argc, char* argv[])
{
	size_t L = 32, U = 1000, S = 1000, seed = 0;
	int imode = 2, filetype = 0;
	std::string gnuplot_args = "''";
	char pottsVersion = 'a';
	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(L, "lattice")["-L"]["--lattice-size"]("Width and height of lattice")
			| Opt(U, "timesteps")["-U"]["--update-steps"]("Number of update steps")
			| Opt(S, "timesteps")["-S"]["--skip-steps"]("Number of updates to skip")
			| Opt(imode, "mode")["--init-mode"]("Initialisation mode")
			| Opt(seed, "seed")["--seed"]("Random seed (0 for unpredictable)")
			| Opt(filetype, "type")["--out-type"]("Output type, 0 - png, 1 - eps")
			| Opt(pottsVersion, "type")["--potts-version"]("Potts version to use. Use b for Ising")
			| Opt(gnuplot_args, "'<args>'")["--gnuplot-args"]("Arguments for gnuplot command");
		auto result = cli.parse(Args(argc, argv));
		if (!result) {
			fmt::print("Error in command line: {}\n", result.errorMessage());
			return EXIT_FAILURE;
		}
		if (!gnuplot_args.empty() && (gnuplot_args.front() != '\'' || gnuplot_args.back() != '\'')) {
			fmt::print("Error in command line: --gnuplot-args must be single-quoted string");
			return EXIT_FAILURE;
		}
		if (showHelp) {
			std::cerr << cli << std::endl;
			return EXIT_SUCCESS;
		}
	}
	const std::string gp_args = gnuplot_args.empty() ? "" : gnuplot_args.substr(1, gnuplot_args.size() - 2);
	const std::string gpcmd = "gnuplot " + gp_args;

	print_debug_info();

	if(pottsVersion != 'a' && pottsVersion != 'b') { ERRPT; fmt::print(stderr, "pottsVersion must be one of 'a' or 'b'"); exit(EXIT_FAILURE); }

	size_t const N = L*L;
	double const inv_L_sq = 1.0 / static_cast<double>(N);
	std::vector<u8> lattice(N);
	std::vector<u8> lattice_buffer(N);
	std::vector<double> energy;
	energy.reserve(U);

	make_dir("frames");

	//double T_min = 0.001;
	//double T_max = 1.5;
	int T_num = 20;
	//double T_step = (T_max - T_min) / (T_num-1.0);

	std::vector<double> avg_energy;
	avg_energy.reserve(T_num);

	std::vector<double> mi;
	mi.reserve(U);

	Hist mi_hist(numStates*numStates, 0);
	std::vector<Hist> prev_mi_hist(10, Hist(numStates * numStates, 1));// (numStates*numStates, 1);
	Hist te_hist(numStates*numStates*numStates, 0);
	Hist gte_hist(static_cast<int64_t>(std::pow(numStates, 6)), 1);
	std::vector<Hist> prev_gte_hist(10, Hist(static_cast<int64_t>(std::pow(numStates, 6)), 1));// std::pow(numStates, 6), 0);
	Hist gte_binary_hist(numStates*numStates * 16, 0);
	Hist gte_reduced_hist(numStates*numStates*numNeighbourhoods, 0);

	Hist mi_uniform(numStates*numStates, 1);
	Hist gte_uniform(static_cast<int64_t>(std::pow(numStates, 6)), 1);

	std::vector<double> kl_mi_uniform;
	std::vector<double> kl_mi_running;
	std::vector<double> kl_gte_uniform;
	std::vector<double> kl_gte_running;

	Hist tmp_mi_hist(numStates*numStates, 0u);
	size_t const kl_measure_length = 10000;
	//for(auto T = T_min; T <= T_max; T += T_step) {
	double T_crit = critical_temp(pottsVersion);
	int ignore_current_energy = 0;
	auto T = 0.5;// T_crit;// 0.3;// 0.375;
	std::mt19937_64 rand_engine(seed ? seed : std::random_device()());
	std::uniform_real_distribution<> rng(0, 1);

	int init = imode;
	if(imode == 3) {
		init = ((T < T_crit) ? 2 : 0);
	}
	initialise(rand_engine, lattice, L, init);
	auto tp_table = tp_table_init(T, pottsVersion);

	FILE* gpipe = open_gnuplot(gpcmd.c_str());

	// Run skip steps
	for(size_t s = 0; s < S; ++s) {
		lattice_buffer = lattice;
		update_glauber(rand_engine, rng, lattice, L, N, tp_table, ignore_current_energy);

		//all_histograms(lattice, lattice_buffer, L, mi_hist, te_hist, gte_hist, gte_binary_hist, gte_reduced_hist, false, rand_engine);

		if((s % kl_measure_length) == 0 && s != 0) {
			kl_mi_uniform.emplace_back(kl_divergence(mi_uniform, mi_hist));
			kl_gte_uniform.emplace_back(kl_divergence(gte_uniform, gte_hist));

			kl_mi_running.emplace_back(kl_divergence(prev_mi_hist[0], mi_hist));
			kl_gte_running.emplace_back(kl_divergence(prev_gte_hist[0], gte_hist));

			for(int64_t i = 0; i < 9; ++i) {
				prev_mi_hist[i] = prev_mi_hist[i + 1];
				prev_gte_hist[i] = prev_gte_hist[i + 1];
			}
			prev_mi_hist[9] = mi_hist;
			prev_gte_hist[9] = gte_hist;

			clear_and_resize_vec(mi_hist);
			clear_and_resize_vec(te_hist);
			clear_and_resize_vec(gte_hist);
			clear_and_resize_vec(gte_binary_hist);
		}
	}
		
	energy.emplace_back(calculate_total_energy(lattice, L, pottsVersion) * inv_L_sq);

	init_frame(gpipe, 2, U+1, L, seed, filetype);
	draw_energy(gpipe, U, energy, mi);
	draw_lattice(gpipe, L, lattice);
	finish_frame(gpipe);

	// Run update steps
	energy.clear();
	for(size_t u = 0; u < U; ++u) {
		update_glauber(rand_engine, rng, lattice, L, N, tp_table, ignore_current_energy);
		energy.emplace_back(calculate_total_energy(lattice, L, pottsVersion) * inv_L_sq);

		init_frame(gpipe, 2, u, L, seed, filetype);
		draw_energy(gpipe, U, energy, mi);
		draw_lattice(gpipe, L, lattice);
		finish_frame(gpipe);
	}

	pclose(gpipe);
	return EXIT_SUCCESS;
}