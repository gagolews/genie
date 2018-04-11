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

#ifndef __HCLUST2_NNBASED_SINGLE_H
#define __HCLUST2_NNBASED_SINGLE_H



// ************************************************************************


#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <deque>
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

   HClustStats stats;

#ifdef _OPENMP
   omp_lock_t pqwritelock;
#endif

   DisjointSets ds;
   bool prefetch;

   virtual void getNearestNeighborsFromMinRadius(size_t index, size_t clusterIndex, double minR, NNHeap& nnheap) = 0;
   void getNearestNeighbors(std::priority_queue<HeapHierarchicalItem> & pq, size_t index);

   void computePrefetch(std::priority_queue<HeapHierarchicalItem> & pq);
   void computeMerge(std::priority_queue<HeapHierarchicalItem> & pq, HClustResult& res);


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
