/* ************************************************************************* *
 *   This file is part of the `DataStructures` package.                      *
 *                                                                           *
 *   Copyright 2015 Maciej Bartoszuk, Anna Cena, Marek Gagolewski,           *
 *                                                                           *
 *   'DataStructures' is free software: you can redistribute it and/or       *
 *   modify it under the terms of the GNU Lesser General Public License      *
 *   as published by the Free Software Foundation, either version 3          *
 *   of the License, or (at your option) any later version.                  *
 *                                                                           *
 *   'DataStructures' is distributed in the hope that it will be useful,     *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU Lesser General Public License for more details.                     *
 *                                                                           *
 *   You should have received a copy of the GNU Lesser General Public        *
 *   License along with 'DataStructures'.                                    *
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
// #include <deque>
// #include <exception>
// #include <string>
// #include <boost/property_map/property_map.hpp>
// #include <boost/tuple/tuple_comparison.hpp>
// #include <algorithm>

#include "hclust2_common.h"
#include "hclust2_merge.h"
#include "disjoint_sets.h"
#include "hclust2_result.h"

namespace DataStructures
{



class HClustNNbasedSingle
{
protected:

   HClustOptions opts;
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

   void computePrefetch(std::priority_queue<HeapHierarchicalItem>& pq);
   void computeMerge(std::priority_queue<HeapHierarchicalItem>& pq, HClustResult& res);


public:

   HClustNNbasedSingle(Distance* dist, RObject control);
   virtual ~HClustNNbasedSingle();

   virtual void print() { Rcout << "this print method is a stub" << std::endl; }

   HClustResult compute();

   inline const HClustStats& getStats()     { return stats; }
   inline const HClustOptions& getOptions() { return opts; }

}; // class

} // namespace DataStructures


#endif
