#ifndef DOS_TESTS
#define DOS_TESTS

#include <catch2/catch.hpp>
#include "../include/model.hpp"
#include "../include/stats.hpp"
#include "../include/defs.hpp"
#include <random>
#include <vector>


double sum(std::vector<double>& joint) {
	return std::accumulate(joint.begin(), joint.end(), 0.);
}

constexpr auto num = numStates * numStates*numSiteEnergy;

TEST_CASE("Normalisation of joint hist", "[NORM]") {
	Hist hist;
	std::vector<double> joint;

	SECTION("Zeroes") {
		hist.resize(64, 0);
		joint.resize(hist.size());
		normalise_hist(hist, joint);
		//Shouldn't do anything
	}

	SECTION("Ones") {
		hist.resize(64, 1);
		joint.resize(hist.size());
		normalise_hist(hist, joint);
		REQUIRE(joint[0] == Approx(1.0 / 64));
		REQUIRE(sum(joint) == Approx(1.0));
	}

	SECTION("Arbitrary 1") {
		hist.resize(64, 1);
		joint.resize(hist.size());
		for(int i = 0; i < 64; ++i) {
			hist[i] = i;
		}
		normalise_hist(hist, joint);
		REQUIRE(sum(joint) == Approx(1.0));
	}

	
	SECTION("Arbitrary 2") {
		
		hist.resize(num, 1);
		joint.resize(hist.size());
		for(int i = 0; i < num; ++i) {
			hist[i] = i;
		}
		normalise_hist(hist, joint);
		REQUIRE(sum(joint) == Approx(1.0));
	}

	SECTION("Random") {
		std::mt19937_64 rand_engine(std::random_device{}());
		std::uniform_int_distribution<> rng(0);
		hist.resize(num);
		joint.resize(hist.size());

		for(auto& h : hist) { h = static_cast<size_t>(rng(rand_engine)); }
		normalise_hist(hist, joint);
		REQUIRE(sum(joint) == Approx(1.0));
	}
}

TEST_CASE("GTE Calculation", "[GTE]") {
	Hist hist;

	SECTION("Zeroes") {
		hist.resize(num, 0);
		//Shouldn't do anything
		REQUIRE(calc_GTE_reduced_from_hist(hist) == Approx(0.0));
	}

	SECTION("Ones") {
		hist.resize(num, 1);
		REQUIRE(calc_GTE_reduced_from_hist(hist) == Approx(0.0));
	}
}

#endif // !DOS_TESTS
