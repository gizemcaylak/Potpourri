/*
* Copyright (C) 2018 Serhan Yýlmaz
*
* This file is part of SPADIS
*
* SPADIS is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SPADIS is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Optimizer.h"
#include <algorithm>

spadis::Optimizer::Optimizer(std::vector<double> scores, Graph graph) : scores_(scores), graph_(graph)
{
	std::vector<size_t> sorted;
	sortedScoreReverseIndices_.resize(scores.size());
	sorted.resize(scores.size());
	std::iota(sorted.begin(), sorted.end(), 0);
	std::sort(sorted.begin(), sorted.end(),
		[&scores](size_t i1, size_t i2) {return scores[i1] > scores[i2]; });
	for (size_t i = 0; i < sorted.size(); i++) {
		sortedScoreReverseIndices_.at(sorted[i]) = i;
	}
	sortedScoreIndices_ = std::move(sorted);
	initialize();
}

void spadis::Optimizer::initialize()
{
	measurements_.BetaMax = 0;
	measurements_.BetaMin = INFINITY;
	flagsDijkstra_.resize(scores_.size());
	handlesDijkstra_.resize(scores_.size());
	for (unsigned int i = 0; i < scores_.size(); i++) {
		flagsDijkstra_.at(i) = Flag::WHITE;
	}
}

void spadis::Optimizer::select(const Options options)
{
	solutions_.clear();
	size_t N = scores_.size();
	std::vector<double> penaltySums;
	std::vector<size_t> betaMaxIndices;
	std::vector<bool> betaMaxFlags;
	if (options.isInMeasurementMode()) {
		betaMaxFlags.resize(scores_.size());
		penaltySums.resize(scores_.size());
	}
	std::vector<OptimizationHandle> optimizationHandles;
	optimizationHandles.resize(scores_.size());
	std::vector<std::vector<Neighbor>> neighborSets;
	neighborSets.reserve(scores_.size());
	for (unsigned int i = 0; i < scores_.size(); i++) {
		neighborSets.push_back(std::vector<Neighbor>());
	}
	double scoreMax = scores_[sortedScoreIndices_[0]];
	for (unsigned int betaIndex = 0; betaIndex < options.getBetaSize(); betaIndex++) {
		size_t lastScoreIndex = 0;
		double lastScore = scoreMax;
		betaMaxIndices.clear();
		optimizationHandles.clear();
		OptimizationQueue optimizationQueue;
		Solution solution;
		const double beta = options.getBeta(betaIndex);
		double betaConstant = 1.0 / (2 * options.getK());
		double betaPrime = beta * betaConstant;
		if (beta == BETA_INFINITE) {
			betaPrime = betaConstant;
		}
		solution.indicators.reserve(N);
		for (size_t i = 0; i < N; i++) {
			solution.indicators.push_back(false);
			OptimizerValue value(scores_.at(i));
			optimizationHandles.push_back(optimizationQueue.push(OptimizerHeapData(value, i)));
		}
		if (options.isInMeasurementMode()) {
			for (size_t i = 0; i < N; i++) {
				betaMaxFlags.at(i) = false;
				penaltySums.at(i) = 0;
			}
		}
		int nSelection = 0;
		OptimizerValue Ftotal;
		while (nSelection < options.getK()) {
			auto a = optimizationQueue.top();
			OptimizerValue Fmax = a.key;
			unsigned int maxFIndex = a.index;
			optimizationQueue.pop();
			auto handlex = optimizationHandles.at(maxFIndex);
			(*handlex).key.infinite -= 1;
			optimizationQueue.update(handlex);
			double scoreCurrent = scores_.at(maxFIndex);
			double penaltyCurrent;
			if (options.isInBetaMeasurementMode()) {
				penaltyCurrent = penaltySums.at(maxFIndex);
				size_t sort_index = sortedScoreReverseIndices_.at(maxFIndex);
				if (beta == BETA_INFINITE) {
					for (auto i = lastScoreIndex + 1; i < sort_index; i++) {
						size_t index = sortedScoreIndices_.at(i);
						if (penaltySums.at(index) > penaltyCurrent
							&& !betaMaxFlags[index]) {
							betaMaxIndices.push_back(index);
							betaMaxFlags[index] = true;
						}
					}
					//if (betaMaxIndices.size() > N) {
					//	DebugUnit::printErrorMessage("Error in betaMaxIndices!");
					//}
					lastScore = scoreCurrent;
					lastScoreIndex = sort_index;
					double betaMax = measurements_.BetaMax;
					for(size_t i = 0; i < betaMaxIndices.size(); i++){
						size_t index = betaMaxIndices[i];
						double score = scores_[index];
						double penalty = penaltySums.at(index);
						if (penaltyCurrent >= penalty)
							continue;
						if (solution.indicators.at(index))
							continue;
						double temp = (score - scoreCurrent) / 
							(betaConstant * (penalty - penaltyCurrent));
						betaMax = std::max(betaMax, temp);
					}
					measurements_.BetaMax = betaMax;
				} else {
					double betaMin = measurements_.BetaMin;
					bool stop = false;
					for (unsigned int i = options.getK(); i < N && !stop; i++) {
						size_t index = sortedScoreIndices_[i];
						double score = scores_[index];
						double penalty = penaltySums.at(index);
						double betaMinPotential = (scoreCurrent - score)
							/ (betaConstant * (penaltyCurrent));
						if (betaMinPotential != 0 && betaMinPotential >= betaMin) {
							stop = true;
						}
						if (penaltyCurrent <= penalty)
							continue;
						double temp = (scoreCurrent - score)
							/ (betaConstant * (penaltyCurrent - penalty));
						betaMin = temp != 0 ? std::min(betaMin, temp) : betaMin;
					}
					measurements_.BetaMin = betaMin;
				}
			}
			if (beta != BETA_INFINITE) {
				Ftotal.real += beta + Fmax.real;
			} else {
				Ftotal.real += Fmax.real;
				Ftotal.infinite += 1 + Fmax.infinite;
			}
			solution.indicators.at(maxFIndex) = true;
			if (neighborSets.at(maxFIndex).empty()) {
				neighborSets.at(maxFIndex) = findNeighbors(maxFIndex, options);
			}
			for (unsigned int j = 0; j < neighborSets.at(maxFIndex).size(); j++) {
				auto pair = neighborSets.at(maxFIndex).at(j);
				double Kvalue = 2 * (1 - pair.second / options.getDistanceParameter());
				if (options.isInBetaMeasurementMode()) {
					penaltySums.at(pair.first) += Kvalue;
					if (beta == BETA_INFINITE) {
						if (scores_.at(pair.first) >= scoreCurrent
							&& penaltySums.at(pair.first) > penaltyCurrent
							&& !betaMaxFlags[pair.first]
							&& !solution.indicators.at(pair.first)) {
							betaMaxFlags[pair.first] = true;
							betaMaxIndices.push_back(pair.first);
						}
					}
				}
				if(!solution.indicators.at(pair.first)) {
					auto handle = optimizationHandles.at(pair.first);
					if (beta != BETA_INFINITE) {
						(*handle).key.real -= betaPrime * Kvalue;
					} else {
						(*handle).key.infinite -= betaPrime * Kvalue;
					}
					optimizationQueue.update(handle);
				}
			}
			nSelection++;
		}
		solution.optimizerValue = Ftotal;
		solutions_.push_back(std::move(solution));
	}
}

std::vector<spadis::Solution> spadis::Optimizer::getSolutions() const
{
	return solutions_;
}

std::pair<double, double> spadis::Optimizer::measureBetaRange(int k, double D)
{
	Options opts;
	opts.setK(k);
	opts.setDistanceParameter(D);
	opts.setBetaMeasurementMode(true);
	opts.addBeta(0);
	opts.addBeta(spadis::BETA_INFINITE);
	select(opts);
	return std::pair<double, double>(measurements_.BetaMin, measurements_.BetaMax);
}

double spadis::Optimizer::measureDeltaRange(int k)
{
	Options opts;
	opts.setK(k);
	opts.setDeltaMeasurementMode(true);
	double Dmin = INFINITY;
	for (size_t i = 0; i < k; i++) {
		size_t index = sortedScoreIndices_.at(i);
		findNeighbors(index, opts, Dmin);
	}
	return Dmin;
}

std::vector<spadis::Optimizer::Neighbor> spadis::Optimizer::findNeighbors(unsigned int sourceIndex, const Options options)
{
	double D = options.getDistanceParameter();
	return std::move(findNeighbors(sourceIndex, options, D));
}

std::vector<spadis::Optimizer::Neighbor> spadis::Optimizer::findNeighbors(unsigned int sourceIndex, const Options options, double& D)
{
	std::vector<Neighbor> out;
	if (!graph_.isWeighted()) {		// Breath-First Search
		std::queue<unsigned int> q;
		std::queue<unsigned int> q_next;
		q.push(sourceIndex);
		int d = 0;
		flagsDijkstra_.at(sourceIndex) = Flag::BLACK;
		flagModificationList_.push_back(sourceIndex);
		while (!q.empty()) {
			if (d >= D) {
				break;
			}
			while (!q.empty()) {
				auto index = q.front();
				if (options.isInDeltaMeasurementMode()) {
					size_t sort_index = sortedScoreReverseIndices_.at(index);
					if (index != sourceIndex && sort_index < options.getK()) {
						D = d;
						q = std::queue<unsigned int>();
						break;
					}
				}
				q.pop();
				out.push_back(Neighbor(index, d));
				const Node node = graph_.getNode(index);
				for (unsigned int i = 0; i < node.size(); i++) {
					unsigned int ix = node.getEdge(i);
					if (flagsDijkstra_.at(ix) == Flag::WHITE) {
						flagsDijkstra_.at(ix) = Flag::BLACK;
						flagModificationList_.push_back(ix);
						q_next.push(ix);
					}
				}
			}
			std::swap(q, q_next);
			d = d + 1;
		}
	} else {						// Dijkstra's Algorithm
		PQDijkstra q;
		handlesDijkstra_.at(sourceIndex) = q.push(DijkstraHeapData(0, sourceIndex));
		flagsDijkstra_.at(sourceIndex) = Flag::BLACK;
		flagModificationList_.push_back(sourceIndex);
		while (!q.empty()) {
			auto a = q.top();
			q.pop();
			if (a.key >= D) {
				break;
			}
			if (options.isInDeltaMeasurementMode()) {
				size_t sort_index = sortedScoreReverseIndices_.at(a.index);
				if (a.index != sourceIndex && sort_index < options.getK()) {
					D = a.key;
					break;
				}
			}
			out.push_back(Neighbor(a.index, a.key));
			flagsDijkstra_.at(a.index) = Flag::BLACK;
			const Node node = graph_.getNode(a.index);
			for (unsigned int i = 0; i < node.size(); i++) {
				unsigned int ix = node.getEdge(i);
				double w = node.getEdgeWeight(i);
				if (flagsDijkstra_.at(ix) == Flag::WHITE) {
					flagsDijkstra_.at(ix) = Flag::GRAY;
					flagModificationList_.push_back(ix);
					handlesDijkstra_.at(ix) = q.push(DijkstraHeapData(a.key + w, ix));
				}
				else if (flagsDijkstra_.at(ix) == Flag::GRAY) {
					auto h = handlesDijkstra_.at(ix);
					if ((*h).key > a.key + w) {
						(*h).key = a.key + w;
						q.increase(h);
					}
				}
			}
		}
	}
	for (int i = 0; i < flagModificationList_.size(); i++) {
		unsigned int index = flagModificationList_.at(i);
		flagsDijkstra_.at(index) = Flag::WHITE;
	}
	flagModificationList_.clear();
	return std::move(out);
}