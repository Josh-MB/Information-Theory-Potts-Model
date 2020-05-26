#ifndef VP_TREE_H
#define VP_TREE_H
// A VP-Tree implementation, by Steve Hanov. (steve.hanov@gmail.com)
// http://stevehanov.ca/blog/index.php?id=130
// This implementation has been modified to use iterative approaches rather than recursive
// Additionally, we preallocate node storage rather than one at a time, since we intend to be
// using a large amount of nodes. Individual allocation results in a significant amount of overhead
// and reduced performance.
// Released to the Public Domain
// Based on "Data Structures and Algorithms for Nearest Neighbor Search" by Peter N. Yianilos
#include "defs.hpp"

#include <fmt/format.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <queue>
#include <limits>
#include <stack>
#include <cmath>
#include <random>
#include <memory>

template<typename T, double(*distance)(const T&, const T&, const double wrapSize)>
class VpTree
{
public:
	VpTree(double const L = TWOPI, std::mt19937_64*const rng = nullptr) : _items(0), L(L), rng(rng), _root(0), rootIdx(-1), nodes(0), nodeCount(0), totalNodes(0) 
	{
		if (!this->rng) {
			std::random_device rd;
			default_rng = std::unique_ptr<std::mt19937_64>(new std::mt19937_64(rd()));
			this->rng = default_rng.get();
		}
	}

	~VpTree() {
		delete[] nodes;
	}

	//Something is wrong in here
	//crashes when shrinking a tree to use with smaller data
	void create(std::vector<T>& items) {
		_items = &items;
		if((int)items.size() > totalNodes) {
			delete[] nodes;
			nodes = 0;
			// fmt::print("Resize tree to {} (from {})\n", items.size(), totalNodes);
			totalNodes = (int)items.size();
			//This approach might run into issues later with HUGE datasets and fragmented memory
			//May have to create a chunked allocator that can allocated large blocks that chain together
			//to fit in fragmented memory more easily.
			nodes = new Node[items.size()];
		}
		nodeCount = 0;
		rootIdx = buildFromPoints(0, (int)items.size());
		_root = &(nodes[rootIdx]);
	}

	/**
	 * Returns the k nearest neighbours to target point and their distances.
	 * This function should be thread safe
	 */
	void knn_search(const T& target, int k, std::vector<T>* results, std::vector<int>* indices,
				std::vector<double>* distances) const
	{
		std::priority_queue<HeapItem> heap;

		double tau = std::numeric_limits<int>::max();
		knnSearchImpl(rootIdx, target, k, heap, tau);

		results->clear(); distances->clear();

		while(!heap.empty()) {
			results->push_back((*_items)[heap.top().index]);
			indices->push_back(heap.top().index);
			distances->push_back(heap.top().dist);
			heap.pop();
		}

		std::reverse(results->begin(), results->end());
		std::reverse(indices->begin(), indices->end());
		std::reverse(distances->begin(), distances->end());
	}
	
	/**
	 * Gets the kth nearest neighbour to target point and its distance
	 * This function should be thread safe
	 */
	void find_kth_neighbour(const T& target, int k, T& kth_neighbour, double& dist) const
	{
		std::priority_queue<HeapItem> heap;

		double tau = std::numeric_limits<int>::max();
		knnSearchImpl(rootIdx, target, k, heap, tau);

		kth_neighbour = (*_items)[heap.top().index];
		dist = heap.top().dist;
	}

	/**
	 * Returns points within radius of target point
	 * This function should be thread safe
	 */
	int fr_search(const T& target, double maxdist, std::vector<T>* results = nullptr) const
	{
		return frSearchImpl(rootIdx, target, maxdist, results);
	}

	int getNodeCount() const { return nodeCount; };
private:
	std::vector<T>* _items;
	const double L;
	std::unique_ptr<std::mt19937_64> default_rng; // Only used if rng isn't provided
	std::mt19937_64* rng;
	std::uniform_real_distribution<double> dis;

	struct Node
	{
		int index;
		double threshold;
		int left;
		int right;

		Node() :
			index(0), threshold(0.), left(-1), right(-1) {}

	}*_root;
	int rootIdx;

	Node* nodes;
	int nodeCount;
	int totalNodes;

	struct HeapItem {
		HeapItem(int index, double dist) :
			index(index), dist(dist) {}
		int index;
		double dist;
		bool operator<(const HeapItem& o) const {
			return dist < o.dist;
		}
	};

	struct DistanceComparator
	{
		const T& item;
		const double L1;
		DistanceComparator(const T& item, const double L) : item(item), L1(L) {}
		bool operator()(const T& a, const T& b) {
			return distance(item, a, L1) < distance(item, b, L1);
		}
	};

	int buildFromPoints(int lower, int upper)
	{
		if(upper == lower) {
			return -1;
		}

		if(nodeCount >= totalNodes) {
			fmt::print("Not enough nodes created, {}/{}\n", nodeCount, totalNodes);
			return -1;
		}
		int nodeId = nodeCount;
		Node* node = &(nodes[nodeCount++]);
		node->index = lower;

		if(upper - lower > 1) {

			// choose an arbitrary point and move it to the start
			double m = dis(*rng);
			int i = (int)(m * (upper - lower - 1)) + lower;
			//int i = (int)(0.5 * (upper - lower - 1)) + lower;
			std::swap((*_items)[lower], (*_items)[i]);

			int median = (upper + lower) / 2;

			// partitian around the median distance
			std::nth_element(
				_items->begin() + lower + 1,
				_items->begin() + median,
				_items->begin() + upper,
				DistanceComparator((*_items)[lower], L));

			// what was the median?
			node->threshold = distance((*_items)[lower], (*_items)[median], L);

			node->index = lower;
			node->left = buildFromPoints(lower + 1, median);
			node->right = buildFromPoints(median, upper);
		}

		return nodeId;
	}

	void knnSearchImpl(int nodeId, const T& target, int k,
				std::priority_queue<HeapItem>& heap, double& tau) const
	{
		if(nodeId == -1) return;
		Node* node = &nodes[nodeId];

		double dist = distance((*_items)[node->index], target, L);

		if(dist < tau) {
			if((int)heap.size() == k) heap.pop();
			heap.push(HeapItem(node->index, dist));
			if((int)heap.size() == k) tau = heap.top().dist;
		}

		if(node->left == -1 && node->right == -1) {
			return;
		}

		if(dist < node->threshold) {
			if(dist - tau <= node->threshold) {
				knnSearchImpl(node->left, target, k, heap, tau);
			}

			if(dist + tau >= node->threshold) {
				knnSearchImpl(node->right, target, k, heap, tau);
			}

		}
		else {
			if(dist + tau >= node->threshold) {
				knnSearchImpl(node->right, target, k, heap, tau);
			}

			if(dist - tau <= node->threshold) {
				knnSearchImpl(node->left, target, k, heap, tau);
			}
		}
	}

	int frSearchImpl(int nodeId, const T& target, const double maxdist, std::vector<T>* results) const
	{
		if(nodeId == -1) return 0;
		//Node* node = &(nodes[nodeId]);

		std::stack<int> workingSet;
		int currNode = -1;
		workingSet.push(nodeId);
		int count = 0;

		while(workingSet.empty() == false || currNode != -1) {
			if(currNode != -1) {
				Node* currNodePtr = &(nodes[currNode]);
				double dist = distance((*_items)[currNodePtr->index], target, L);

				if(dist < maxdist) {
					++count;
					if(results) results->push_back((*_items)[currNodePtr->index]);
				}
				if(currNodePtr->left == -1 && currNodePtr->right == -1) {
					currNode = -1;
					currNodePtr = NULL;
					continue;
				}

				if(dist < currNodePtr->threshold) {
					if((dist - maxdist <= currNodePtr->threshold))
						workingSet.push(currNodePtr->left);
					if((dist + maxdist >= currNodePtr->threshold))
						workingSet.push(currNodePtr->right);
				}
				else {
					if((dist + maxdist >= currNodePtr->threshold))
						workingSet.push(currNodePtr->right);
					if((dist - maxdist <= currNodePtr->threshold))
						workingSet.push(currNodePtr->left);
				}

			}
			if(workingSet.empty() == false) {
				currNode = workingSet.top();
				workingSet.pop();
			}
		}

		return count;
	}
};

#endif //VP_TREE_H
