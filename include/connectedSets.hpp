#ifndef CONNECTED_SETS_HPP_
#define CONNECTED_SETS_HPP_

#include "defs.hpp"

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <queue>
#include <functional>

// Calculates the connected sets in a particle cloud
struct ConnectedSets
{
	//Set of connected sets
	//Inner vector defines a connected set of particles (id's are stored)
	//Outer vector defines the set of connected sets.
	//sum(inner.size()) = N
	//Might turn this into list<vector<int>>. Could be more efficient. Means don't have to preallocate
	//everything, and can deal with holes much easier
	std::vector<std::vector<int> > sets;
	//Map mapping each particle to the set it is part of
	std::vector<int> particleMap;
	std::priority_queue<int, std::vector<int>, std::greater<int> > holes;
	int nextSet; //Last active set

	// Initialise to handle N particles
	ConnectedSets(size_t N) {
		sets.resize(N);
		for (size_t i = 0; i < N; ++i) {
			sets[i].reserve(N);
		}
		particleMap.resize(N, -1);
		nextSet = 0;
	}

	ConnectedSets(const ConnectedSets& rhs) {
		sets.resize(rhs.sets.size());
		for (size_t i = 0; i < sets.size(); ++i) {
			sets[i].assign(rhs.sets[i].begin(), rhs.sets[i].end());
		}
		particleMap.assign(rhs.particleMap.begin(), rhs.particleMap.end());
		nextSet = rhs.nextSet;
		holes = rhs.holes;
	}

	ConnectedSets& operator=(const ConnectedSets& rhs) {
		sets.resize(rhs.sets.size());
		for (size_t i = 0; i < sets.size(); ++i) {
			sets[i].assign(rhs.sets[i].begin(), rhs.sets[i].end());
		}
		particleMap.assign(rhs.particleMap.begin(), rhs.particleMap.end());
		nextSet = rhs.nextSet;
		holes = rhs.holes;
		return *this;
	}

	void clear() {
		for (size_t i = 0; i < sets.size(); ++i) {
			sets[i].clear();
		}

		for (size_t i = 0; i < particleMap.size(); ++i) {
			particleMap[i] = -1;
		}
		nextSet = 0;
	}

	//Connects two elements ids. If both are already in the map,
	//will assign all from upper set to lower set. If one is
	//not in map, it will be added to other's set. Otherwise
	//new set is created with both elements in it
	void connect(int idx_i, int idx_j) {
		//Make new set with both particles
		if (particleMap[idx_i] == -1 && particleMap[idx_j] == -1) {
			int newSet = -1;
			if (!holes.empty()) {
				newSet = holes.top();
				holes.pop();
			}
			else {
				newSet = nextSet++;
			}
			particleMap[idx_i] = particleMap[idx_j] = newSet;
			sets[newSet].push_back(idx_i);
			sets[newSet].push_back(idx_j);
		}
		//Merge sets into lower
		else if (particleMap[idx_i] != -1 && particleMap[idx_j] != -1) {
			int lowerSet = std::min(particleMap[idx_i], particleMap[idx_j]);
			int upperSet = std::max(particleMap[idx_i], particleMap[idx_j]);

			if (lowerSet != upperSet) {
				if (upperSet == nextSet)
					nextSet--;
				else
					holes.push(upperSet);

				for (size_t i = 0; i < sets[upperSet].size(); ++i) {
					sets[lowerSet].push_back(sets[upperSet][i]);
					particleMap[sets[upperSet][i]] = lowerSet;
				}
				sets[upperSet].clear();
			}
		}
		//Add idx_i into idx_j's set
		else if (particleMap[idx_i] == -1)
		{
			particleMap[idx_i] = particleMap[idx_j];
			sets[particleMap[idx_j]].push_back(idx_i);
		}
		//Add idx_j into idx_i's set
		else if (particleMap[idx_j] == -1)
		{
			particleMap[idx_j] = particleMap[idx_i];
			sets[particleMap[idx_i]].push_back(idx_j);
		}
	}

	//Add an element to the map. Generally only do this if it is not
	//connected to any other elements. If it is connected, just use
	//the connect function (even if neither are already in the map)
	void add(int idx_i) {
		if (particleMap[idx_i] == -1) {
			particleMap[idx_i] = nextSet;
			sets[nextSet].push_back(idx_i);
			nextSet++;
		}

	}

	//Rearranges sets to leave no holes
	void squash()
	{
		while (!holes.empty())
		{
			int nextHole = holes.top();
			holes.pop();

			int lastSet = nextSet - 1;
			for (; lastSet >= 0 && sets[lastSet].empty(); lastSet--) {
			}
			if (lastSet < nextHole) {
				//all holes have been filled. Only holes exist are past end of last full set
				while (!holes.empty())
					holes.pop();
				for (; nextSet > 0 && sets[nextSet - 1].empty(); nextSet--) {}
				break;
			}

			for (size_t i = 0; i < sets[lastSet].size(); ++i) {
				sets[nextHole].push_back(sets[lastSet][i]);
				particleMap[sets[lastSet][i]] = nextHole;
			}
			sets[lastSet].clear();
			nextSet--;
		}
	}
};

#endif