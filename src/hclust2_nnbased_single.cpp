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


#include "hclust2_nnbased_single.h"

using namespace Rcpp;
using namespace std;
using namespace boost;
using namespace DataStructures;


// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustNNbasedSingle::HClustNNbasedSingle(Distance* dist, RObject control) :
      opts(control),
      _n(dist->getObjectCount()),
      _distance(dist),
      _indices(dist->getObjectCount()),
      // _indicesinv(dist->getObjectCount()),
      neighborsCount(dist->getObjectCount(), 0),
      minRadiuses(dist->getObjectCount(), -INFINITY),
      // maxRadiuses(vector<double>(dist->getObjectCount(), INFINITY)),
      shouldFind(dist->getObjectCount(), true),
      nearestNeighbors(dist->getObjectCount()),
   #ifdef GENERATE_STATS
      stats(),
   #endif
      ds(dist->getObjectCount())
{
   // starting _indices: random permutation of {0,1,...,_n-1}
   for (size_t i=0;i<_n;i++)
      _indices[i] = i;
   for (size_t i=_n-1; i>= 1; i--)
      swap(_indices[i], _indices[(size_t)(unif_rand()*(i+1))]);
}


HClustNNbasedSingle::~HClustNNbasedSingle() {

}



HeapNeighborItem HClustNNbasedSingle::getNearestNeighbor(size_t index, double distMax)
{
   size_t clusterIndex = ds.find_set(index);
   if (shouldFind[index] && nearestNeighbors[index].empty())
   {
      if (minRadiuses[index] > distMax) {
         return HeapNeighborItem(SIZE_MAX, minRadiuses[index]);
      }

#ifdef GENERATE_STATS
#ifdef _OPENMP
#pragma omp atomic
#endif
      ++stats.nnCals;
#endif
      NNHeap nnheap((prefetch)?opts.maxNNPrefetch:opts.maxNNMerge);
      getNearestNeighborsFromMinRadius(index, clusterIndex, minRadiuses[index], nnheap);
      nnheap.fill(nearestNeighbors[index]);

      size_t newNeighborsCount = nearestNeighbors[index].size();

      neighborsCount[index] += newNeighborsCount;
      if (neighborsCount[index] > _n - index || newNeighborsCount == 0)
         shouldFind[index] = false;

      if (newNeighborsCount > 0)
         minRadiuses[index] = nearestNeighbors[index].back().dist;
   }

   if (!nearestNeighbors[index].empty())
   {
      while (!nearestNeighbors[index].empty()) {
#ifdef GENERATE_STATS
#ifdef _OPENMP
#pragma omp atomic
#endif
         ++stats.nnCount;
#endif
         auto res = nearestNeighbors[index].front();
         nearestNeighbors[index].pop_front();
         if (clusterIndex != ds.find_set(res.index))
            return res;
      }
      return HeapNeighborItem(SIZE_MAX, minRadiuses[index]);
   }
   else
   {
      return HeapNeighborItem(SIZE_MAX, INFINITY);
   }
}



void HClustNNbasedSingle::computePrefetch(std::priority_queue<HeapHierarchicalItem>& pq)
{
   // INIT: Pre-fetch a few nearest neighbors for each point
   MESSAGE_2("[%010.3f] prefetching NNs\n", clock()/(float)CLOCKS_PER_SEC);

#ifdef _OPENMP
   omp_set_dynamic(0); /* the runtime will not dynamically adjust the number of threads */
   omp_lock_t writelock;
   omp_init_lock(&writelock);
   #pragma omp parallel for schedule(dynamic)
#endif
   for (size_t i=0; i<_n; i++)
   {
#ifndef _OPENMP
      Rcpp::checkUserInterrupt(); // may throw an exception, fast op, not thread safe
#endif
      HeapNeighborItem hi=getNearestNeighbor(i);
      if (hi.index != SIZE_MAX)
      {
#if !defined(_OPENMP)
         MESSAGE_7("\r             prefetch NN: %d/%d", i, _n-1);
#endif
#ifdef _OPENMP
         omp_set_lock(&writelock);
#endif
         pq.push(HeapHierarchicalItem(i, hi.index, hi.dist));
#ifdef _OPENMP
         omp_unset_lock(&writelock);
#endif
      }
   }
#ifdef _OPENMP
   omp_destroy_lock(&writelock);
#endif
   MESSAGE_7("\r             prefetch NN: %d/%d\n", _n-1, _n-1);
}


void HClustNNbasedSingle::computeMerge(
      std::priority_queue<HeapHierarchicalItem>& pq,
      HClustResult& res)
{
   MESSAGE_2("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);

   size_t i = 0;
   while (true)
   {
      HeapHierarchicalItem hhi = pq.top();
      pq.pop();

      if (hhi.index2 == SIZE_MAX) {
         HeapNeighborItem hi=getNearestNeighbor(hhi.index1, INFINITY);
         if (isfinite(hi.dist))
            pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));
         continue;
      }

      size_t s1 = ds.find_set(hhi.index1);
      size_t s2 = ds.find_set(hhi.index2);
      if (s1 != s2)
      {
         Rcpp::checkUserInterrupt(); // may throw an exception, fast op

         res.link(_indices[hhi.index1], _indices[hhi.index2], hhi.dist);
         ds.link(s1, s2);

         ++i;
         if (i == _n-1) break; /* avoid computing unnecessary nn */
      }
      MESSAGE_7("\r             %d / %d", i+1, _n);

      STOPIFNOT(hhi.index1 < hhi.index2);
      HeapNeighborItem hi=getNearestNeighbor(hhi.index1, pq.top().dist);
      STOPIFNOT(hhi.index1 < hi.index);
      if (isfinite(hi.dist))
         pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));
   }
   MESSAGE_7("\r             %d / %d\n", _n, _n);
   Rcpp::checkUserInterrupt();
}


HClustResult HClustNNbasedSingle::compute()
{
   priority_queue<HeapHierarchicalItem> pq;
   HClustResult res(_n, _distance);

#if VERBOSE >= 5
   _distance->getStats().print();
#endif

   prefetch = true;
   computePrefetch(pq);
   prefetch = false;

#if VERBOSE >= 5
   _distance->getStats().print();
#endif

   computeMerge(pq, res);

   return res;
}
