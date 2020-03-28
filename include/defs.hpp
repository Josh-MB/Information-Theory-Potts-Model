#pragma once
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <iostream>
#include <numeric>

using u8 = uint8_t;

//Set to 5 since we are only using the 5-state model
u8 constexpr numStates = 5;

int constexpr CURRENT_DOS_FILE_VER = 3;
int constexpr CURRENT_GTE_FILE_VER = 4;

int constexpr numNeighbourhoods = 9; // Used for aggregating neighbours in GTE reduced
int constexpr numSiteEnergy = 5; // Number of possible energy values for a given site (-4, -3, -2, -1, 0)
constexpr int neighbour_sites = 4;

using Hist = std::vector<size_t>;

//Helper macros for wrapping lattice coordinates
inline size_t wrap_minus(size_t u, size_t L1) {return (u == 0) ? L1 : u - 1; }; // -1 with lattice wrap
inline size_t wrap_plus(size_t u, size_t L1) {return (u == L1) ? 0 : u + 1; }; // +1 with lattice wrap

//Helper functions for accumulating and clearing histograms
template <class T>
inline void clear_and_resize_vec(std::vector<T> & v) {auto sz = v.size(); v.clear(); v.resize(sz, 0); };

inline void acc_hist(Hist& source, Hist& dest) {std::transform(source.begin(), source.end(), dest.begin(), dest.begin(), std::plus<size_t>()); };

inline void acc_thread_hists(std::vector<Hist>& threadHists, Hist& dest) {
	for (auto& v : threadHists) {
		std::transform(v.begin(), v.end(), dest.begin(), dest.begin(), std::plus<size_t>());
		clear_and_resize_vec(v);
	}
}

inline int sgn(int const val) {
	return (0 < val) - (val < 0);
}

inline int make_idx(int const x, const int y, const int L) {
	return y * L + x;
}

// Wraps x in boundaries of L. i.e returns
// 0 if x == L, L-1 if x==-1 or x otherwise
inline int wrap(int const x, int const L) {
	return x + ((x == L)*-L) + ((x == -1)*L);
}

//Works out the wrapped difference
inline int unwrap(int const x, int const L) {
	if(x == L - 1)
		return -1;
	else if(x == 1 - L)
		return 1;
	else
		return x;
}

//Helpers for write data to binary file
#define BIN_WRITE(var, type, num, file) if(fwrite(&var, sizeof(type), num, file) != num) {PEEXIT("write failed");}
#define BIN_WRITE_PTR(var, type, num, file) if(fwrite(var, sizeof(type), num, file) != num) {PEEXIT("write failed");}
#define BIN_READ(var, type, num, file) if(fread(&var, sizeof(type), num, file) != num) {PEEXIT("read failed");}
#define BIN_READ_PTR(var, type, num, file) if(fread(var, sizeof(type), num, file) != num) {PEEXIT("read failed");}
#define CHECK_VAL(var, expected, flag) if(var != expected) { std::cout << "vars not equal " << #var << ":" << var << ", " << #expected << ":" << expected << "\n"; flag = false;}

struct DoSFile {
	int fileVersion{ -1 };
	size_t L{ 0 };
	size_t U{ 0 };
	size_t seed{ 0 };
	u8 numStates{ 0 };
	double factor{ 0.0 };
	size_t iteration{ 0 };
	// logG with all elements where count=0 is removed
	std::vector<double> logG_compressed;
	std::vector<double> logG;
	Hist hist;
	std::vector<size_t> indices;
	std::vector<u8> lattice;
};

void read_dos_file(char const*const dosFile, std::vector<double> &g, std::vector<int> &e, size_t& dosL);
// Read DoS data from file. If rand_engine is not nullptr, will update with the saved state
DoSFile read_dos_file(char const* const dosFile, std::mt19937_64* rand_engine = nullptr);

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(),
		 [&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });

	return idx;
}