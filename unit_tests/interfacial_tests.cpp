#ifndef INTERFACIAL_TESTS
#define INTERFACIAL_TESTS

#include <catch2/catch.hpp>
#include "../include/stats.hpp"
#include "../include/defs.hpp"
#include <vector>

TEST_CASE("Basic Interfacial Length", "[IFACE_LEN]") {
	int L = 8;
	std::vector<u8> lattice(L*L, 0);

	SECTION("No elements") {
		auto lengths = calcInterfaceAllLengths(lattice, L);
		REQUIRE(lengths.empty());
	}

	SECTION("Single element - Middle") {
		lattice[make_idx(1, 1, L)] = 1;
		auto lengths = calcInterfaceAllLengths(lattice, L);
		REQUIRE(lengths.size() == 2);
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
	}

	SECTION("Single element - Left Edge") {
		lattice[make_idx(0, 0, L)] = 1;
		auto lengths = calcInterfaceAllLengths(lattice, L);
		REQUIRE(lengths.size() == 2);
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
	}


	SECTION("Two adjacent elements") {
		lattice[make_idx(0, 0, L)] = 1;
		lattice[make_idx(1, 0, L)] = 1;
		auto lengths = calcInterfaceAllLengths(lattice, L);
		REQUIRE(lengths.size() == 2);
		REQUIRE(lengths[0] == 6);
		REQUIRE(lengths[1] == 6);
	}

	SECTION("Two diagonal elements") {
		lattice[make_idx(0, 0, L)] = 1;
		lattice[make_idx(1, 1, L)] = 1;
		auto lengths = calcInterfaceAllLengths(lattice, L);
		REQUIRE(lengths.size() == 3);
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
		REQUIRE(lengths[2] == 8);
	}


	SECTION("Two Off diagonal elements") {
		lattice[make_idx(1, 0, L)] = 1;
		lattice[make_idx(0, 1, L)] = 1;
		auto lengths = calcInterfaceAllLengths(lattice, L);
		REQUIRE(lengths.size() == 3);
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
		REQUIRE(lengths[2] == 8);
	}

	SECTION("Vertical Strip") {
		for(int y = 0; y < L; ++y) {
			lattice[make_idx(1, y, L)] = 1;
		}
		auto lengths = calcInterfaceAllLengths(lattice, L);
		REQUIRE(lengths.size() == 4);
		REQUIRE(lengths[0] == L);
		REQUIRE(lengths[1] == L);
		REQUIRE(lengths[2] == L);
		REQUIRE(lengths[3] == L);
	}

	SECTION("Diagonal Strip") {
		for(int y = 0; y < L; ++y) {
			lattice[make_idx(y, y, L)] = 1;
		}
		auto lengths = calcInterfaceAllLengths(lattice, L);
		REQUIRE(lengths.size() == (L+2));
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
		REQUIRE(lengths[2] == 4);
		REQUIRE(lengths[3] == 4);
		REQUIRE(lengths[4] == 4);
		REQUIRE(lengths[5] == 4);
		REQUIRE(lengths[6] == 4);
		REQUIRE(lengths[7] == 4);
		REQUIRE(lengths[8] == (L * 2));
		REQUIRE(lengths[9] == (L * 2));
	}

	SECTION("Non-Full Diagonal Strip") {
		for(int y = 0; y < L-1; ++y) {
			lattice[make_idx(y, y, L)] = 1;
		}
		auto lengths = calcInterfaceAllLengths(lattice, L);
		REQUIRE(lengths.size() == L);
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
		REQUIRE(lengths[2] == 4);
		REQUIRE(lengths[3] == 4);
		REQUIRE(lengths[4] == 4);
		REQUIRE(lengths[5] == 4);
		REQUIRE(lengths[6] == 4);
		REQUIRE(lengths[7] == ((L-1) * 4));
	}

	SECTION("Off Diagonal Strip") {
		for(int y = 0; y < L; ++y) {
			lattice[make_idx(L-y-1, y, L)] = 1;
		}
		auto lengths = calcInterfaceAllLengths(lattice, L);
		REQUIRE(lengths.size() == (L+2));
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
		REQUIRE(lengths[2] == 4);
		REQUIRE(lengths[3] == 4);
		REQUIRE(lengths[4] == 4);
		REQUIRE(lengths[5] == 4);
		REQUIRE(lengths[6] == 4);
		REQUIRE(lengths[7] == 4);
		REQUIRE(lengths[8] == (L * 2));
		REQUIRE(lengths[9] == (L * 2));

	}

	SECTION("Donut") {
		lattice[make_idx(0, 0, L)] = 1;
		lattice[make_idx(1, 0, L)] = 1;
		lattice[make_idx(2, 0, L)] = 1;
		lattice[make_idx(0, 1, L)] = 1;
		lattice[make_idx(2, 1, L)] = 1;
		lattice[make_idx(0, 2, L)] = 1;
		lattice[make_idx(1, 2, L)] = 1;
		lattice[make_idx(2, 2, L)] = 1;
		auto lengths = calcInterfaceAllLengths(lattice, L);
		REQUIRE(lengths.size() == 4);
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
		REQUIRE(lengths[2] == 12);
		REQUIRE(lengths[3] == 12);
	}

	SECTION("Not Quite Donut") {
		lattice[make_idx(0, 0, L)] = 1;
		lattice[make_idx(1, 0, L)] = 1;
		lattice[make_idx(2, 0, L)] = 1;
		lattice[make_idx(0, 1, L)] = 1;
		lattice[make_idx(2, 1, L)] = 1;
		lattice[make_idx(0, 2, L)] = 1;
		lattice[make_idx(1, 2, L)] = 1;
		auto lengths = calcInterfaceAllLengths(lattice, L);
		REQUIRE(lengths.size() == 3);
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 12);
		REQUIRE(lengths[2] == 16);
	}
}

TEST_CASE("Alternate Interfacial Length", "[IFACE_LEN_ALT]") {
	int L = 8;
	std::vector<u8> lattice(L*L, 0);

	SECTION("No elements") {
		auto lengths = calcInterfaceAllLengthsAlt(lattice, L);
		REQUIRE(lengths.size() == 1);
		REQUIRE(lengths[0] == 0);
	}

	SECTION("Single element - Middle") {
		lattice[make_idx(1, 1, L)] = 1;
		auto lengths = calcInterfaceAllLengthsAlt(lattice, L);
		REQUIRE(lengths.size() == 2);
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
	}

	SECTION("Single element - Left Edge") {
		lattice[make_idx(0, 0, L)] = 1;
		auto lengths = calcInterfaceAllLengthsAlt(lattice, L);
		REQUIRE(lengths.size() == 2);
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
	}


	SECTION("Two adjacent elements") {
		lattice[make_idx(0, 0, L)] = 1;
		lattice[make_idx(1, 0, L)] = 1;
		auto lengths = calcInterfaceAllLengthsAlt(lattice, L);
		REQUIRE(lengths.size() == 2);
		REQUIRE(lengths[0] == 6);
		REQUIRE(lengths[1] == 6);
	}

	SECTION("Two diagonal elements") {
		lattice[make_idx(0, 0, L)] = 1;
		lattice[make_idx(1, 1, L)] = 1;
		auto lengths = calcInterfaceAllLengthsAlt(lattice, L);
		REQUIRE(lengths.size() == 3);
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
		REQUIRE(lengths[2] == 8);
	}


	SECTION("Two Off diagonal elements") {
		lattice[make_idx(1, 0, L)] = 1;
		lattice[make_idx(0, 1, L)] = 1;
		auto lengths = calcInterfaceAllLengthsAlt(lattice, L);
		REQUIRE(lengths.size() == 3);
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
		REQUIRE(lengths[2] == 8);
	}

	SECTION("Vertical Strip") {
		for(int y = 0; y < L; ++y) {
			lattice[make_idx(1, y, L)] = 1;
		}
		auto lengths = calcInterfaceAllLengthsAlt(lattice, L);
		REQUIRE(lengths.size() == 2);
		REQUIRE(lengths[0] == 2*L);
		REQUIRE(lengths[1] == 2*L);
	}

	SECTION("Diagonal Strip") {
		for(int y = 0; y < L; ++y) {
			lattice[make_idx(y, y, L)] = 1;
		}
		auto lengths = calcInterfaceAllLengthsAlt(lattice, L);
		REQUIRE(lengths.size() == (L + 1));
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
		REQUIRE(lengths[2] == 4);
		REQUIRE(lengths[3] == 4);
		REQUIRE(lengths[4] == 4);
		REQUIRE(lengths[5] == 4);
		REQUIRE(lengths[6] == 4);
		REQUIRE(lengths[7] == 4);
		REQUIRE(lengths[8] == (L * 4));
	}

	SECTION("Non-Full Diagonal Strip") {
		for(int y = 0; y < L - 1; ++y) {
			lattice[make_idx(y, y, L)] = 1;
		}
		auto lengths = calcInterfaceAllLengthsAlt(lattice, L);
		REQUIRE(lengths.size() == L);
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
		REQUIRE(lengths[2] == 4);
		REQUIRE(lengths[3] == 4);
		REQUIRE(lengths[4] == 4);
		REQUIRE(lengths[5] == 4);
		REQUIRE(lengths[6] == 4);
		REQUIRE(lengths[7] == ((L - 1) * 4));
	}

	SECTION("Off Diagonal Strip") {
		for(int y = 0; y < L; ++y) {
			lattice[make_idx(L - y - 1, y, L)] = 1;
		}
		auto lengths = calcInterfaceAllLengthsAlt(lattice, L);
		REQUIRE(lengths.size() == (L + 1));
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 4);
		REQUIRE(lengths[2] == 4);
		REQUIRE(lengths[3] == 4);
		REQUIRE(lengths[4] == 4);
		REQUIRE(lengths[5] == 4);
		REQUIRE(lengths[6] == 4);
		REQUIRE(lengths[7] == 4);
		REQUIRE(lengths[8] == (L * 4));

	}

	SECTION("Donut") {
		lattice[make_idx(0, 0, L)] = 1;
		lattice[make_idx(1, 0, L)] = 1;
		lattice[make_idx(2, 0, L)] = 1;
		lattice[make_idx(0, 1, L)] = 1;
		lattice[make_idx(2, 1, L)] = 1;
		lattice[make_idx(0, 2, L)] = 1;
		lattice[make_idx(1, 2, L)] = 1;
		lattice[make_idx(2, 2, L)] = 1;
		auto lengths = calcInterfaceAllLengthsAlt(lattice, L);
		REQUIRE(lengths.size() == 3);
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 12);
		REQUIRE(lengths[2] == 16);
	}

	SECTION("Not Quite Donut") {
		lattice[make_idx(0, 0, L)] = 1;
		lattice[make_idx(1, 0, L)] = 1;
		lattice[make_idx(2, 0, L)] = 1;
		lattice[make_idx(0, 1, L)] = 1;
		lattice[make_idx(2, 1, L)] = 1;
		lattice[make_idx(0, 2, L)] = 1;
		lattice[make_idx(1, 2, L)] = 1;
		auto lengths = calcInterfaceAllLengthsAlt(lattice, L);
		REQUIRE(lengths.size() == 3);
		std::sort(lengths.begin(), lengths.end());
		REQUIRE(lengths[0] == 4);
		REQUIRE(lengths[1] == 12);
		REQUIRE(lengths[2] == 16);
	}
}


#endif //INTERFACIAL_TESTS