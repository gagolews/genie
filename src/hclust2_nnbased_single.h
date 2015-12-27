/* ************************************************************************* *
 *   This file is part of the `grup` package.                                *
 *                                                                           *
 *   Copyright 2015 Maciej Bartoszuk, Anna Cena, Marek Gagolewski,           *
 *                                                                           *
 *   'grup' is free software: you can redistribute it and/or                 *
 *   modify it under the terms of the GNU Lesser General Public License      *
 *   as published by the Free Software Foundation, either version 3          *
 *   of the License, or (at your option) any later version.                  *
 *                                                                           *
 *   'grup' is distributed in the hope that it will be useful,               *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU Lesser General Public License for more details.                     *
 *                                                                           *
 *   You should have received a copy of the GNU Lesser General Public        *
 *   License along with 'grup'.                                              *
 *   If not, see <http://www.gnu.org/licenses/>.                             *
 * ************************************************************************* */

#ifndef __HCLUST2_NNBASED_SINGLE_H
#define __HCLUST2_NNBASED_SINGLE_H



// ************************************************************************


#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
// #include <fstream>
#include <deque>
// #include <exception>
// #include <string>
// #include <boost/property_map/property_map.hpp>
// #include <boost/tuple/tuple_comparison.hpp>
// #include <algorithm>

#include "hclust2_common.h"
#include "disjoint_sets.h"
#include "hclust2_result.h"

using namespace std;
using namespace Rcpp;

namespace grup
{



class HClustNNbasedSingle
{
protected:

   HClustOptions* opts;
   size_t n;
   Distance* distance;
   std::vector<size_t> indices;

   std::vector<size_t> neighborsCount;
   std::vector<double> minRadiuses;
   std::vector<bool> shouldFind;
   std::vector< deque<HeapNeighborItem> > nearestNeighbors;

   HClustStats stats;

   DisjointSets ds;
   bool prefetch;

   virtual void getNearestNeighborsFromMinRadius(size_t index, size_t clusterIndex, double minR, NNHeap& nnheap) = 0;
   HeapNeighborItem getNearestNeighbor(size_t index, double distMax=INFINITY);

   void computePrefetch(HclustPriorityQueue& pq);
   void computeMerge(HclustPriorityQueue& pq, HClustResult& res);


public:

   HClustNNbasedSingle(Distance* dist, HClustOptions* opts);
   virtual ~HClustNNbasedSingle();

   virtual void print() { Rcout << "this print method is a stub" << std::endl; }

   HClustResult compute(bool lite=false);

   inline const HClustStats& getStats()     { return stats; }
   inline const HClustOptions& getOptions() { return *opts; }

}; // class

} // namespace grup


#endif
