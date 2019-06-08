#include "../include/scopeTimer.hpp"
#include <fmt/format.h>
#include <clara.hpp>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <cstring>

int sim_anim(int argc, char* argv[]);
int sim_glauber(int argc, char* argv[]);
int sim_kl_filtered_gte(int argc, char* argv[]);
int sim_dos(int argc, char* argv[]);
int sim_dos_gte_et(int argc, char* argv[]);
int sim_order(int argc, char* argv[]);
int sim_dos_quantities(int argc, char* argv[]);
int sim_record_data(int argc, char* argv[]);

int main(int argc, char* argv[])
{
	if(!(argc>1)) {
		fmt::print(stderr, "{}: must be at least one program argument\n", argv[0]);
		return EXIT_FAILURE;
	}

	bool showHelp = false;
	std::string command;
	using namespace clara::detail;
	auto cli = Help(showHelp)
		| Arg(command, "submodule")("Which submodule to tun:\n"
			"anim - Perform visualisation\n"
			"glauber - Calculate metrics using Glauber updating\n"
			"glauber_kl - As with 'glauber', but filtering GTE histogram based on KL divergence\n"
			"dos - Calculate the Density of State\n"
			"dos_gte - Calculate GTE(E,T) from DoS to derive GTE from in matlab\n"
			"thermodynamic_quantities - Calculate thermodynamic MI, magnetisation and interface lengths from DoS\n"
			"record_data - Perform Glauber simulation and write all lattice states to file for external processing\n"
			"order - Calculate magnetisation and autocorrelation time from Glauber simulation");
	auto result = cli.parse(Args(2, argv)); // Just read the first two args
	if (!result) {
		fmt::print("Error in command line: {}\n", result.errorMessage());
		return EXIT_FAILURE;
	}
	if (showHelp) {
		std::cerr << cli << std::endl;
		return EXIT_SUCCESS;
	}

	int retCode = EXIT_FAILURE;
	{
		ScopeTimer timer("Potts");
		int argc_1 = argc - 1;	char** argv_1 = argv + 1; // Clara expects just one ignored argument, so move forward by one
		if (command == "anim") retCode = sim_anim(argc_1, argv_1);
		else if (command == "glauber") retCode = sim_glauber(argc_1, argv_1);
		else if (command == "glauber_kl") retCode = sim_kl_filtered_gte(argc_1, argv_1);
		else if (command == "dos") retCode = sim_dos(argc_1, argv_1);
		else if (command == "dos_gte") retCode = sim_dos_gte_et(argc_1, argv_1);
		else if (command == "order") retCode = sim_order(argc_1, argv_1);
		else if (command == "thermodynamic_quantities") retCode = sim_dos_quantities(argc_1, argv_1);
		else if (command == "record_data") retCode = sim_record_data(argc_1, argv_1);
		else {
			fmt::print(stderr, "{}: unknown submodule '{}'\n", argv[0], argv[1]);
			return EXIT_FAILURE;
		}
	}

    return retCode;
}