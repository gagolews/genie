/* ************************************************************************* *
 *   This file is part of the `genie` package for R.                         *
 *                                                                           *
 *   Copyright 2015-2018 Marek Gagolewski, Maciej Bartoszuk, Anna Cena       *
 *                                                                           *
 *   'genie' is free software: you can redistribute it and/or                *
 *   modify it under the terms of the GNU General Public License             *
 *   as published by the Free Software Foundation, either version 3          *
 *   of the License, or (at your option) any later version.                  *
 *                                                                           *
 *   'genie' is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with 'genie'. If not, see <http://www.gnu.org/licenses/>.         *
 * ************************************************************************* */

#ifndef __HCLUST2_NNBASED_GINI_H
#define __HCLUST2_NNBASED_GINI_H



// ************************************************************************


#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <deque>
#include "hclust2_common.h"
#include "disjoint_sets.h"
#include "hclust2_result.h"

namespace grup
{



class HClustNNbasedGini
{
protected:

   HClustOptions* opts;
   size_t n;
   Distance* distance;
   std::vector<size_t> indices;

   std::vector<size_t> neighborsCount;
   std::vector<double> minRadiuses;
   std::vector<bool> shouldFind;
   std::vector< std::deque<HeapNeighborItem> > nearestNeighbors;

   HClustStats stats;

   bool symmetric;

   PhatDisjointSets ds;
   bool prefetch;

   virtual void getNearestNeighborsFromMinRadius(size_t index, size_t clusterIndex, double minR, double& maxR, NNHeap& nnheap) = 0;
   HeapNeighborItem getNearestNeighbor(size_t index, double distMax=INFINITY);

   void prefetchNNsSymmetric();
   void computePrefetch(HclustPriorityQueue& pq);
   void computeMerge(HclustPriorityQueue& pq, HClustResult& res);
   void linkAndRecomputeGini(double& lastGini, size_t s1, size_t s2);

public:

   HClustNNbasedGini(Distance* dist, HClustOptions* opts);
   virtual ~HClustNNbasedGini();

   HClustResult compute();

   inline const HClustStats& getStats()     { return stats; }
   inline const HClustOptions& getOptions() { return *opts; }

}; // class

} // namespace grup


#endif
