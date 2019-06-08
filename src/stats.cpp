#include "../include/stats.hpp"
#include "../include/model.hpp"
#include "../include/utils.hpp"

#include <algorithm>
#include <cassert>

void MI_histogram(std::vector<u8> const & lattice, size_t const L, Hist& hist)
{
	if(hist.size() != numStates*numStates) PEEXIT("MI histogram is not correct size");
	if(lattice.size() != L*L) PEEXIT("Lattice size or L is wrong");

	size_t const L1 = L - 1;
	for(auto y = 0u; y < L; ++y) {
		for(auto x = 0u; x < L; ++x) {
			//Only look up "right" and "down" neighbour,
			//but add both combinations for each, i.e (i,r), (r,i), (i,d), (d,i)
			//Means less lookups
			int64_t s_i = lattice[y*L + x];
			int64_t s_r = lattice[y*L + wrap_plus(x, L1)];
			int64_t s_d = lattice[wrap_plus(y, L1)*L + x];
			++hist[s_i*numStates + s_r];
			++hist[s_i*numStates + s_d];
			++hist[s_r*numStates + s_i];
			++hist[s_d*numStates + s_i];
		}
	}

	/*size_t const L1 = L - 1;
	for(auto y = 0u; y < L; ++y) {
		for(auto x = 0u; x < L; ++x) {
			auto e = -calc_site_energy(lattice, x, y, L);
			auto s = lattice[y*L + x];
			++hist[s*numStates + e];
		}
	}*/
}

void TE_histogram(std::vector<u8> const & lattice, std::vector<u8> const & lattice_old, size_t const L, Hist& hist)
{
	if(hist.size() != numStates*numStates*numStates) PEEXIT("TE histogram is not correct size");
	if(lattice.size() != L*L) PEEXIT("Lattice size or L is wrong");
	if(lattice_old.size() != L*L) PEEXIT("Old lattice size or L is wrong");


	size_t const L1 = L - 1;
	for(auto y = 0u; y < L; ++y) {
		for(auto x = 0u; x < L; ++x) {
			int64_t s_i_new = lattice[y*L + x];
			int64_t s_i = lattice_old[y*L + x];
			int64_t s_u = lattice_old[wrap_minus(y, L1)*L + x];
			int64_t s_d = lattice_old[wrap_plus(y, L1)*L + x];
			int64_t s_l = lattice_old[y*L + wrap_minus(x, L1)];
			int64_t s_r = lattice_old[y*L + wrap_plus(x, L1)];

			size_t s_i_idx = s_i_new + s_i*numStates;
			size_t fact = numStates*numStates;
			++hist[s_i_idx + s_u*fact];
			++hist[s_i_idx + s_d*fact];
			++hist[s_i_idx + s_l*fact];
			++hist[s_i_idx + s_r*fact];
		}
	}
}

void GTE_histogram(std::vector<u8> const & lattice, std::vector<u8> const & lattice_old, size_t const L, Hist& hist)
{
	if(hist.size() != static_cast<unsigned int>(std::pow(numStates, 6))) PEEXIT("GTE histogram is not correct size");
	//static_assert(std::pow(numStates, 6) == 15625);
	if(lattice.size() != L*L) PEEXIT("Lattice size or L is wrong");
	if(lattice_old.size() != L*L) PEEXIT("Old lattice size or L is wrong");


	size_t const L1 = L - 1;
	for(auto y = 0u; y < L; ++y) {
		for(auto x = 0u; x < L; ++x) {
			int64_t s_i_new = lattice[y*L + x];
			int64_t s_i = lattice_old[y*L + x];
			int64_t s_u = lattice_old[wrap_minus(y, L1)*L + x];
			int64_t s_d = lattice_old[wrap_plus(y, L1)*L + x];
			int64_t s_l = lattice_old[y*L + wrap_minus(x, L1)];
			int64_t s_r = lattice_old[y*L + wrap_plus(x, L1)];

			//size_t s_i_idx = s_i_new + numStates * s_i*;
			++hist[s_i_new + numStates * (s_i +
				numStates * (s_u +
				numStates * (s_d +
				numStates * (s_l +
				numStates * s_r))))];
		}
	}
}

void GTE_reduced_histogram(std::vector<u8> const & lattice, std::vector<u8> const & lattice_old, size_t const L, Hist& hist)
{
	// Not actually needed (encoded in next table), but used to remember which
	// element is which
	//int bin_lookup[] = {	0,	/*AAAA*/	1,	/*AAAX*/	2,	/*AAXX*/	3,	/*AXXX*/	4,	/*XXXX*/
	//						5,	/*AAAW*/	6,	/*AAWX*/	7,	/*AWXX*/	8,	/*WXXX*/
	//						9,	/*AAWW*/	10,	/*AWWX*/	11,	/*WWXX*/
	//						12,	/*AWWW*/	13,	/*WWWX*/
	//						14,	/*WWWW*/ };
	// Lookup which gives the bin number for a given (E_X,E_W) pair. -1 values are for X,W combinations that are
	// not possible, i.e E_X=E_W=4.
	// Note if X=W, then setE_W to 0 and increment (E_X, 0) and (0, E_X) instead
	//constexpr int xw_lookup[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, -1, 9, 10, 11, -1, -1, 12, 13, -1, -1, -1, 14, -1, -1, -1, -1 };
	//std::uniform_int_distribution<> toggle(0, 1);

	if(hist.size() != numStates*numStates*numSiteEnergy) PEEXIT("Reduced GTE histogram is not correct size");
	if(lattice.size() != L * L) PEEXIT("Lattice size or L is wrong");
	if(lattice_old.size() != L * L) PEEXIT("Old lattice size or L is wrong");

	size_t const L1 = L - 1;
	for(auto y = 0u; y < L; ++y) {
		for(auto x = 0u; x < L; ++x) {
			int64_t s_i_new = lattice[y*L + x];
			int64_t s_i = lattice_old[y*L + x];
			int64_t s_u = lattice_old[wrap_minus(y, L1)*L + x];
			int64_t s_d = lattice_old[wrap_plus(y, L1)*L + x];
			int64_t s_l = lattice_old[y*L + wrap_minus(x, L1)];
			int64_t s_r = lattice_old[y*L + wrap_plus(x, L1)];

			int e_x = (s_i == s_u) + (s_i == s_d) + (s_i == s_l) + (s_i == s_r);
			++hist[s_i_new + numStates * (s_i + numStates * static_cast<int64_t>(e_x))];
		}
	}
}

void all_histograms(std::vector<u8> const & lattice, std::vector<u8> const & lattice_old, size_t const L, Hist& mi, Hist& te, Hist& gte, Hist& binary_gte, Hist& reduced_gte, bool sortStates,
					std::mt19937_64 & engine)
{
	if(mi.size() != numStates * numStates) PEEXIT("MI histogram is not correct size");
	if(te.size() != numStates * numStates*numStates) PEEXIT("TE histogram is not correct size");
	if(binary_gte.size() != numStates * numStates * 16) PEEXIT("Binary GTE histogram is not correct size");
	if(reduced_gte.size() != numStates * numStates*numSiteEnergy) PEEXIT("Reduced GTE histogram is not correct size");
	if(gte.size() != static_cast<unsigned int>(std::pow(numStates, 6))) PEEXIT("GTE histogram is not correct size");
	if(lattice.size() != L * L) PEEXIT("Lattice size or L is wrong");
	if(lattice_old.size() != L * L) PEEXIT("Old lattice size or L is wrong");
	
	//If T>=T_c, and E<-1.25*N
	//do sort order 
	Hist sort_order(numStates, 0);
	std::iota(sort_order.begin(), sort_order.end(), 0);
	if(sortStates) {
		Hist counts(numStates, 0);
		calc_magnetisation_and_counts(lattice, L, counts);

		std::uniform_int_distribution<> flip(0, 1);
		std::uniform_int_distribution<> state(0, numStates - 1);

		std::vector<double> freq(numStates, 0);
		std::transform(counts.begin(), counts.end(), freq.begin(), [L](size_t c) {return static_cast<double>(c) / (L*L); });

		// Randomly swap elements if within 5% (of total) of each other
		for(int i = 0; i < numStates - 1; ++i) {
			int s1 = state(engine);
			int s2 = state(engine);
			if(std::abs(freq[s1] - freq[s2]) < 0.05) {
				if(flip(engine) == 1) {
					std::swap(freq[s1], freq[s2]);
					std::swap(counts[s1], counts[s2]);
				}
			}
		}

		auto temp_sort_order = sort_indexes(counts);
				
		for(int i = 0; i < numStates; ++i) {
			sort_order[temp_sort_order[i]] = i;
		}
	}
	size_t const L1 = L - 1;
	for(auto y = 0u; y < L; ++y) {
		for(auto x = 0u; x < L; ++x) {
			int64_t s_i_new = sort_order[lattice[y*L + x]];
			//u8 s_r_new = lattice[y*L + wrap_plus(x, L1)];
			//u8 s_d_new = lattice[wrap_plus(y, L1)*L + x];
			int64_t s_i = sort_order[lattice_old[y*L + x]];
			int64_t s_u = sort_order[lattice_old[wrap_minus(y, L1)*L + x]];
			int64_t s_d = sort_order[lattice_old[wrap_plus(y, L1)*L + x]];
			int64_t s_l = sort_order[lattice_old[y*L + wrap_minus(x, L1)]];
			int64_t s_r = sort_order[lattice_old[y*L + wrap_plus(x, L1)]];

			int energy = (s_i == s_u) + (s_i == s_d) + (s_i == s_r) + (s_i == s_l);
			
			/*++mi[s_i_new*numStates + s_r_new];
			++mi[s_i_new*numStates + s_d_new];
			++mi[s_r_new*numStates + s_i_new];
			++mi[s_d_new*numStates + s_i_new];*/
			++mi[s_i*numStates + s_r];
			++mi[s_i*numStates + s_d];
			++mi[s_r*numStates + s_i];
			++mi[s_d*numStates + s_i];

			/*size_t s_i_idx = s_i_new + s_i*numStates;
			size_t fact = numStates*numStates;
			++te[s_i_idx + s_u*fact];
			++te[s_i_idx + s_d*fact];
			++te[s_i_idx + s_l*fact];
			++te[s_i_idx + s_r*fact];*/

			++gte[s_i_new + numStates * (s_i +
				  numStates * (s_u +
				  numStates * (s_d +
				  numStates * (s_l +
			      numStates * s_r))))];

			/*++binary_gte[s_i_new + numStates * (s_i + 
				numStates * ((s_i==s_u) + 
				2 * ((s_i==s_d) +
				2 * ((s_i==s_l) +
				2 * ((s_i==s_r))))))];*/

			++reduced_gte[s_i_new + numStates * (s_i + numStates * static_cast<int64_t>(energy))];
		}
	}
}

double calc_MI(Hist const & hist)
{
	std::vector<double> joint_dist(hist.size(), 0.);
	normalise_hist(hist, joint_dist);
	auto count = std::accumulate(hist.begin(), hist.end(), size_t{});
	if(count == 0u) return 0;
	auto sanity_check = std::accumulate(joint_dist.begin(), joint_dist.end(), 0.);
	if(std::abs(sanity_check - 1) > 0.0001) { printf("MI Joint dist doesn't sum to 1\n"); }

	std::vector<double> dist_x(numStates, 0.); //the x marginal (i.e. summed over y)
	std::vector<double> dist_y(numStates, 0.); //the y marginal (i.e. summed over x)

	for(uint64_t y = 0u; y < numStates; ++y) {
		for(uint64_t x = 0u; x < numStates; ++x) {
			dist_y[y] += joint_dist[y*numStates + x];
			dist_x[x] += joint_dist[y*numStates + x];
		}
	}

	return entropy(dist_y) + entropy(dist_x) - entropy(joint_dist);
}

double calc_GTE_from_hist(Hist const & hist)
{
	std::vector<double> joint_dist(hist.size(), 0.);
	//printf("Normalising\n");
	normalise_hist(hist, joint_dist);
	//printf("Constructing marginals\n");

	std::vector<double> xw_dist(numStates*numStates, 0.);
	std::vector<double> xy_dist(static_cast<int64_t>(std::pow(numStates, 5)), 0.);
	std::vector<double> x_dist(numStates, 0.);
	//printf("Counting marginals\n");

	for(uint64_t w = 0u; w < numStates; ++w) {
		for(uint64_t x = 0u; x < numStates; ++x) {
			for(uint64_t u = 0u; u < numStates; ++u) {
				for(uint64_t d = 0u; d < numStates; ++d) {
					for(uint64_t l = 0u; l < numStates; ++l) {
						for(uint64_t r = 0u; r < numStates; ++r) {
							auto idx = w + numStates * (x +
									   numStates * (u +
									   numStates * (d +
									   numStates * (l +
									   numStates * r))));
							auto xyIdx = (x +
										  numStates * (u +
										  numStates * (d +
										  numStates * (l +
										  numStates * r))));

							xw_dist[w + x * numStates] += joint_dist[idx];
							xy_dist[xyIdx] += joint_dist[idx];
							x_dist[x] += joint_dist[idx];
						}
					}
				}
			}
		}
	}
	//printf("Calculating final GTE\n");
	return entropy(xw_dist) + entropy(xy_dist) - entropy(x_dist) - entropy(joint_dist);
	//++gte[s_i_new + numStates * (s_i +
	//	    numStates * (s_u +
	//	    numStates * (s_d +
	//	    numStates * (s_l +
	//	    numStates * s_r))))];
}

double calc_GTE_reduced_from_hist(Hist const & hist)
{
	if(hist.size() != numStates * numStates * numSiteEnergy) PEEXIT("Reduced GTE histogram is not correct size");

	std::vector<double> joint_dist(hist.size(), 0.);
	//printf("Normalising\n");
	normalise_hist(hist, joint_dist);
	auto count = std::accumulate(hist.begin(), hist.end(), size_t{});
	if(count == 0u) return 0; // No elements were in hist, so GTE is 0
	auto sanity_check = std::accumulate(joint_dist.begin(), joint_dist.end(), 0.);
	if(std::abs(sanity_check - 1) > 0.0001) { printf("GTE Joint dist doesn't sum to 1\n"); }
	//printf("Constructing marginals\n");

	std::vector<double> xw_dist(numStates*numStates, 0.);
	std::vector<double> xy_dist(numStates*numSiteEnergy, 0.);
	std::vector<double> x_dist(numStates, 0.);
	//printf("Counting marginals\n");

	for(int64_t w = 0u; w < numStates; ++w) {
		for(int64_t x = 0u; x < numStates; ++x) {
			for(int64_t y = 0u; y < numSiteEnergy; ++y) {
				auto idx = w + numStates * (x + numStates * y);
				auto xyIdx = x + numStates * y;

				xw_dist[w + x * numStates] += joint_dist[idx];
				xy_dist[xyIdx] += joint_dist[idx];
				x_dist[x] += joint_dist[idx];
			}
		}
	}
	
	//printf("Calculating final GTE\n");
	return entropy(xw_dist) + entropy(xy_dist) - entropy(x_dist) - entropy(joint_dist);
	//++gte[s_i_new + numStates * (s_i +
	//	    numStates * (s_u +
	//	    numStates * (s_d +
	//	    numStates * (s_l +
	//	    numStates * s_r))))];
}

double calc_GTE_binary_from_hist(Hist const & hist)
{
	if(hist.size() != numStates * numStates * 16) PEEXIT("Binary GTE histogram is not correct size");

	std::vector<double> joint_dist(hist.size(), 0.);
	//printf("Normalising\n");
	normalise_hist(hist, joint_dist);
	auto count = std::accumulate(hist.begin(), hist.end(), size_t{});
	if(count == 0u) return 0; // No elements were in hist, so GTE is 0
	auto sanity_check = std::accumulate(joint_dist.begin(), joint_dist.end(), 0.);
	if(std::abs(sanity_check - 1) > 0.0001) { printf("GTE Binary Joint dist doesn't sum to 1\n"); }
	//printf("Constructing marginals\n");

	std::vector<double> xw_dist(numStates*numStates, 0.);
	std::vector<double> xy_dist(numStates*16, 0.);
	std::vector<double> x_dist(numStates, 0.);
	//printf("Counting marginals\n");

	for(int64_t w = 0u; w < numStates; ++w) {
		for(int64_t x = 0u; x < numStates; ++x) {
			for(int64_t u = 0u; u < 2; ++u) {
				for(int64_t d = 0u; d < 2; ++d) {
					for(int64_t l = 0u; l < 2; ++l) {
						for(int64_t r = 0u; r < 2; ++r) {
							int64_t idx = w + numStates * (x +
													numStates * (u +
													2 * (d +
													2 * (l +
													2 * r))));
							auto xyIdx = (x +
										  numStates * (u +
													   2 * (d +
																	2 * (l +
																				 2 * r))));

							xw_dist[w + x * numStates] += joint_dist[idx];
							xy_dist[xyIdx] += joint_dist[idx];
							x_dist[x] += joint_dist[idx];
						}
					}
				}
			}
		}
	}

	//printf("Calculating final GTE\n");
	return entropy(xw_dist) + entropy(xy_dist) - entropy(x_dist) - entropy(joint_dist);
	//++gte[s_i_new + numStates * (s_i +
	//	    numStates * (s_u +
	//	    numStates * (s_d +
	//	    numStates * (s_l +
	//	    numStates * s_r))))];
}

double entropy(std::vector<double> const & dist)
{
	auto ret = 0.;
	for(auto i : dist)
	{
		ret += ((i == 0) ? 0 : i* std::log(i));
	}
	return -ret;
}

void normalise_hist(Hist const& hist, std::vector<double> & normalised_hist)
{
	if(hist.size() != normalised_hist.size()) PEEXIT("Hist sizes in normalise_hist are different");
	auto count = std::accumulate(hist.begin(), hist.end(), 0.);
	if(count == 0) return;
	auto total_inv = 1.0 / count;
	std::transform(hist.begin(), hist.end(), normalised_hist.begin(), [total_inv](size_t l) {return l * total_inv; });

}

double kl_divergence(Hist const& baseHist, Hist const&  targetHist)
{
	if(baseHist.size() != targetHist.size()) PEEXIT("Hist sizes in KL divergence calculation are different");

	auto baseTotal = static_cast<double>(std::accumulate(baseHist.begin(), baseHist.end(), size_t{}));
	auto targetTotal = static_cast<double>(std::accumulate(targetHist.begin(), targetHist.end(), size_t{}));

	size_t size = baseHist.size();

	double invBaseTotal = 1. / baseTotal;
	double invTargetTotal = 1. / targetTotal;

	double kl = 0.;
	for(auto i = 0u; i < size; ++i)
	{
		double px = baseHist[i] * invBaseTotal;
		size_t t = (targetHist[i] != 0) ? targetHist[i] : 1;
		double logTerm = 0;
		if(baseHist[i] != 0)
			logTerm = std::log2(px / (t * invTargetTotal));
		kl += px * logTerm;
	}

	return kl;
}

int scaleEnergyValue(int const energy, size_t const L) {
	if(L < 32) PEEXIT("Scaling only works with L>=32)");
	if(energy == 0) return 0;
	return (energy - 1) / static_cast<int>((L / 32) * (L / 32)) + 1;
}

void turnRight(int& fx, int& fy) {
	int tmp = sgn(fy); fy = sgn(fx); fx = -tmp;
}

void turnLeft(int& fx, int& fy) {
	int tmp = sgn(fy); fy = -sgn(fx); fx = tmp;
}

enum EdgeMasks
{
	EM_RIGHT = 1,
	EM_DOWN = 2,
	EM_LEFT = 4,
	EM_UP = 8,
};

void applyMask(int px, int py, int nx, int ny, int const L, std::vector<u8> &edgeMask) {
	/*int const L1 = L - 1;
	auto tmpx = std::minmax(px, nx);
	int x = ((tmpx.first == 0 && tmpx.second == L1) ? tmpx.second : tmpx.first);
	auto tmpy = std::minmax(py, ny);
	int y = ((tmpy.first == 0 && tmpy.second == L1) ? tmpy.second : tmpy.first);*/
	auto i = make_idx(px, py, L); //i is wrong on border
								  // px,py should be adjacent to nx, ny, thus one of these abs calls will be zero
								  // the sgn check is to handle boundary wrapping
								  //edgeMask[i] |= EM_RIGHT * std::abs(sgn(px - nx)) + EM_DOWN * std::abs(sgn(py - ny));
	if(px == nx) {
		if(py == ny + 1 || (py == 0 && ny == (L - 1))) {
			edgeMask[i] |= EM_UP;
			return;
		}
		if(py == ny - 1 || (py == (L - 1) && ny == 0)) {
			edgeMask[i] |= EM_DOWN;
			return;
		}
		printf("Shouldn't get here (up/down): p=(%d, %d), n=(%d, %d)\n", px, py, nx, ny);
	}
	else if(py == ny) {
		if(px == nx + 1 || (px == 0 && nx == (L - 1))) {
			edgeMask[i] |= EM_LEFT;
			return;
		}
		if(px == nx - 1 || (px == (L - 1) && nx == 0)) {
			edgeMask[i] |= EM_RIGHT;
			return;
		}
		printf("Shouldn't get here (left/right): p=(%d, %d), n=(%d, %d)\n", px, py, nx, ny);
	}
	printf("Shouldn't get here (outer): p=(%d, %d), n=(%d, %d)\n", px, py, nx, ny);
}
int walkInterface(std::vector<u8> const & lattice, int px, int py, int nx, int ny, u8 s_i, int const L,
				  std::vector<u8> &cellUsed, std::vector<u8> &edgeMasks) {
	int length = 0;
	int const origPx = px, origPy = py, origNx = nx, origNy = ny;
	do {
		//Apply mask, increment length, mark the cell as used, then crawl one step
		//if(px == origPx && py == origPy)
			applyMask(px, py, nx, ny, L, edgeMasks);
		int pidx = make_idx(px, py, L);
		cellUsed[pidx] = 1;
		++length;

		// If pn = p-n, then forward vec is (pn.y, -pn.x)
		// use sgn to normalise
		int fx = unwrap(py - ny, L);
		int fy = unwrap(nx - px, L);

		//Work out front right and front left cells
		int frx = wrap(nx + fx - sgn(fy), L);
		int fry = wrap(ny + fy + sgn(fx), L);
		int frIdx = make_idx(frx, fry, L);

		int flx = wrap(px + fx + sgn(fy), L);
		int fly = wrap(py + fy - sgn(fx), L);
		int flIdx = make_idx(flx, fly, L);

		// If front left cell is same, turn left, regardless of what right cell is
		// note neighbour shouldn't be equal to current cell
		if(lattice[frIdx] != s_i) {
			// Note turn right, means p is stationary
			nx = frx;
			ny = fry;
			continue;
		}

		// If front right is same (given that front left is not),
		// then go straight, otherwise turn right
		if(lattice[flIdx] != s_i)
		{
			// Both elements move forward one spot, which is just
			// the forward right/left cells
			px = frx;
			py = fry;
			nx = flx;
			ny = fly;
			continue;
		}

		// n is stationary, p changes
		px = flx;
		py = fly;
	} while(!(px == origPx && py == origPy && nx == origNx && ny == origNy));

	return length;
}

/**
* Calculates the interface length of every cluster in the lattice
* Does this by iterating through every (interface) edge and walking along until
* cycle is found, accumulating length as it goes
*/
std::vector<int> calcInterfaceAllLengths(std::vector<u8> const & lattice, int const L)
{
	std::vector<int> lengths;
	int const L1 = L - 1;
	// Bit is set if the right neighbour edge has been walked, or
	// down neighbour edge has been walked (To check other edges,
	// check that neighbour's right/down edge)
	std::vector<u8> edgeMasks(lattice.size(), 0);

	// Marks clusters already checked
	std::vector<u8> cellUsed(lattice.size(), 0);

	for(int y = 0; y < L; ++y) {
		for(int x = 0; x < L; ++x) {
			auto const i = make_idx(x, y, L);

			//if(cellUsed[i] == 0) {
				auto const s_i = lattice[i];

				int lx = static_cast<int>(wrap_minus(x, L1));
				int rx = static_cast<int>(wrap_plus(x, L1));
				int uy = static_cast<int>(wrap_minus(y, L1));
				int dy = static_cast<int>(wrap_plus(y, L1));
				int lidx = make_idx(lx, y, L);
				int ridx = make_idx(rx, y, L);
				int uidx = make_idx(x, uy, L);
				int didx = make_idx(x, dy, L);

				// Look for first neighbour that is different. Check all neighbours as
				// could have two interfaces (i.e donut)
				//should I be checking edgeMasks?
				if(!(edgeMasks[i] & EM_RIGHT) && lattice[ridx] != s_i) {
					lengths.emplace_back(walkInterface(lattice, x, y, rx, y, s_i, L, cellUsed, edgeMasks));
				}
				if(!(edgeMasks[i] & EM_DOWN) && lattice[didx] != s_i) {
					lengths.emplace_back(walkInterface(lattice, x, y, x, dy, s_i, L, cellUsed, edgeMasks));
				}

				// check left and up
				if(!(edgeMasks[i] & EM_LEFT) && lattice[lidx] != s_i) {
					lengths.emplace_back(walkInterface(lattice, x, y, lx, y, s_i, L, cellUsed, edgeMasks));
				}
				if(!(edgeMasks[i] & EM_UP) && lattice[uidx] != s_i) {
					lengths.emplace_back(walkInterface(lattice, x, y, x, uy, s_i, L, cellUsed, edgeMasks));
				}
			//}
		}
	}

	return lengths;
}

std::vector<int> calcInterfaceAllLengthsAlt(std::vector<u8> const& lattice, int const L)
{
	std::vector<int> mask(lattice.size(), -1);

	int numClusters = 0;
	int const L1 = L - 1;
	using Coord = std::pair<int, int>;
	std::vector<Coord> indicesToCheck;
	indicesToCheck.reserve(lattice.size());
	for(int y = 0; y < L; ++y) {
		for(int x = 0; x < L; ++x) {
			auto i = make_idx(x, y, L);
			u8 const s = lattice[i];
			// Skip if already processed
			if(mask[i] != -1)
				continue;

			mask[i] = numClusters;
			indicesToCheck.push_back({ x,y });
			while(!indicesToCheck.empty())
			{
				auto c = indicesToCheck.back();
				i = make_idx(c.first, c.second, L);
				indicesToCheck.pop_back();
				// Calculate 4 neigbour coords
				int x2[] = { static_cast<int>(wrap_plus(c.first, L1)), static_cast<int>(wrap_minus(c.first, L1)),
					c.first, c.first };
				int y2[] = { c.second, c.second, static_cast<int>(wrap_plus(c.second, L1)), 
					static_cast<int>(wrap_minus(c.second, L1)) };
				// Add them to the set (and mask) if same state and not processed
				for(int n = 0; n < 4; ++n) {
					int other = make_idx(x2[n], y2[n], L);
					if(mask[other] == -1 && lattice[other] == s) {
						indicesToCheck.push_back({ x2[n], y2[n] });
						mask[other] = numClusters;
					}
				}

			}

			indicesToCheck.clear();
			++numClusters;
		}
	}

	// Mask now contains all clusters of states. A cluster is all sites s_j which are
	// connected by a path of bonds to s_i, where all states are s_j = s_i
	// Now need to count the perimeter of each cluster (number of edges whose neighbours
	// aren't both s_i)
	std::vector<int> lengths(numClusters, 0);
	for(int y = 0; y < L; ++y) {
		for(int x = 0; x < L; ++x) {
			// Neighbour coords
			int x2[] = { static_cast<int>(wrap_plus(x, L1)), static_cast<int>(wrap_minus(x, L1)), x, x };
			int y2[] = { y, y, static_cast<int>(wrap_plus(y, L1)), static_cast<int>(wrap_minus(y, L1)) };
			auto i = make_idx(x, y, L);
			auto cluster = mask[i];
			// For each neighbour, increment length if not part of same cluster
			for(int n = 0; n < 4; ++n) {
				int other = make_idx(x2[n], y2[n], L);
				lengths[cluster] += (mask[other] != mask[i]);
			}
		}
	}

	return lengths;
}
