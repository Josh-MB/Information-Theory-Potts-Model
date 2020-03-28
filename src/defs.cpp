#include "../include/defs.hpp"
#include "../include/utils.hpp"

#include <stdio.h>
#include <sstream>

void read_dos_file(char const*const dosFile, std::vector<double>& g, std::vector<int>& e, size_t &L)
{
	FILE* fp = fopen(dosFile, "rb");
	if(fp == nullptr) PEEXIT("Failed to open dos file %s", dosFile);

	if(check_magic_header(fp)) {
		printf("dos file magic is good\n");
		size_t L1, U, seed, len, it, sz;
		u8 n;
		int tnum;
		double tmin, tmax, factor;
		int fileVer;

		BIN_READ(fileVer, int, 1, fp);
		BIN_READ(L1, size_t, 1, fp);
		BIN_READ(U, size_t, 1, fp);
		BIN_READ(seed, size_t, 1, fp);
		if(fileVer > 1) {
			BIN_READ(len, size_t, 1, fp);
			printf("allocing rngBuf, len: %zu\n", len);
			char* rngBuf = new char[len];
			BIN_READ_PTR(rngBuf, char, len, fp);
			delete[] rngBuf;
		}
		BIN_READ(n, u8, 1, fp);
		printf("DoS file is using q=%d, L=%zu\n", n, L1);
		L = L1;
		if (fileVer < 3)
		{
			BIN_READ(tnum, int, 1, fp);
			BIN_READ(tmin, double, 1, fp);
			BIN_READ(tmax, double, 1, fp);
		}
		if(fileVer > 1) {
			BIN_READ(factor, double, 1, fp);
			BIN_READ(it, size_t, 1, fp);
		}
		BIN_READ(sz, size_t, 1, fp);
		g.resize(sz);
		e.resize(sz);
		BIN_READ_PTR(&(g[0]), double, sz, fp);
		printf("allocing finalHist\n");
		size_t* finalHist = new size_t[sz];
		BIN_READ_PTR(finalHist, size_t, sz, fp);
		for (size_t i = 0; i < sz; ++i)
		{
			if (finalHist[i] == 0)
				g[i] = -1;
		}
		delete[] finalHist;
		Hist e1(sz, 0);
		BIN_READ_PTR(&(e1[0]), size_t, sz, fp);
		for(int i = 0; i < sz; ++i) {
			e[i] = - static_cast<int>(e1[i]) / 2;
		}
	}
	
}

DoSFile read_dos_file(char const* const dosFile, std::mt19937_64* rand_engine)
{
	DoSFile ret;
	FILE* binFile = fopen(dosFile, "rb");
	if (binFile == nullptr) PEEXIT("Failed to open previous DoS File: %s", dosFile);

	if (!check_magic_header(binFile)) PEEXIT("Bad magic");

	BIN_READ(ret.fileVersion, int, 1, binFile);
	BIN_READ(ret.L, size_t, 1, binFile);
	BIN_READ(ret.U, size_t, 1, binFile);
	BIN_READ(ret.seed, size_t, 1, binFile);
	if (ret.fileVersion > 1) {
		if (rand_engine) {
			rand_engine->seed(ret.seed);
		}
		size_t len;
		BIN_READ(len, size_t, 1, binFile);
		char* buf1 = new char[len + 1];
		BIN_READ_PTR(buf1, char, len, binFile);
		std::istringstream sbuf(buf1);
		std::istream is(sbuf.rdbuf());
		if (rand_engine) {
			is >> *rand_engine;
		}
		delete[] buf1;
	}
	BIN_READ(ret.numStates, u8, 1, binFile);
	if (ret.fileVersion > 1) {
		BIN_READ(ret.factor, double, 1, binFile);
		BIN_READ(ret.iteration, size_t, 1, binFile);
	}

	bool valid = true;
	CHECK_VAL(ret.numStates, numStates, valid);
	if (!valid)
		PEEXIT("Input file is not valid");

	size_t vec_size = ret.L * ret.L * 2 + 1;
	ret.logG.resize(vec_size, 0);
	ret.logG_compressed.resize(vec_size, 0);
	ret.hist.resize(vec_size, 0);
	ret.indices.resize(vec_size, -1);
	size_t L2 = ret.L * ret.L;
	ret.lattice.resize(L2);

	size_t sz = 0, fullSz = vec_size;
	BIN_READ(sz, size_t, 1, binFile);
	if(sz != fullSz) PEEXIT("Data wrong size in binary file")
	BIN_READ_PTR(ret.logG_compressed.data(), double, sz, binFile);
	BIN_READ_PTR(ret.hist.data(), size_t, sz, binFile);
	BIN_READ_PTR(ret.indices.data(), size_t, sz, binFile);
	BIN_READ_PTR(ret.lattice.data(), u8, L2, binFile);

	for (size_t j = 0; j < fullSz; ++j)
	{
		if (ret.indices[j] != -1) {

			if (ret.indices[j] / 2 > fullSz) PEEXIT("Indices are wrong: %d", ret.indices[j] / 2);
			ret.logG[ret.indices[j] / 2] = ret.logG_compressed[j];
		}
	}
	return ret;
}
