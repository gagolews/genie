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


#include "hclust2_nnbased_single.h"

using namespace grup;


// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustNNbasedSingle::HClustNNbasedSingle(Distance* dist, HClustOptions* opts) :
      opts(opts),
      n(dist->getObjectCount()),
      distance(dist),
      indices(dist->getObjectCount()),
      neighborsCount(dist->getObjectCount(), 0),
      minRadiuses(dist->getObjectCount(), -INFINITY),
      shouldFind(dist->getObjectCount(), true),
   #ifdef GENERATE_STATS
      stats(),
   #endif
      ds(dist->getObjectCount())
{
   // starting indices: random permutation of {0,1,...,_n-1}
   for (size_t i=0;i<n;i++)
      indices[i] = i;
   for (size_t i=n-1; i>= 1; i--)
      swap(indices[i], indices[(size_t)(unif_rand()*(i+1))]);

#ifdef _OPENMP
   omp_init_lock(&pqwritelock);
#endif
}


HClustNNbasedSingle::~HClustNNbasedSingle() {
#ifdef _OPENMP
   omp_destroy_lock(&pqwritelock);
#endif
}



void HClustNNbasedSingle::getNearestNeighbors(
   std::priority_queue< HeapHierarchicalItem > & pq,
   size_t index)
{
   if (!shouldFind[index])
      return;

   size_t clusterIndex = ds.find_set(index);
#ifdef GENERATE_STATS
#ifdef _OPENMP
#pragma omp atomic
#endif
   ++stats.nnCals;
#endif
   NNHeap nnheap;
   getNearestNeighborsFromMinRadius(index, clusterIndex, minRadiuses[index], nnheap);
   size_t newNeighborsCount = 0.0;

#ifdef _OPENMP
   omp_set_lock(&pqwritelock);
#endif
   while (!nnheap.empty()) {
      if (isfinite(nnheap.top().dist) && nnheap.top().index != SIZE_MAX) {
         ++newNeighborsCount;
         pq.push(HeapHierarchicalItem(index, nnheap.top().index, nnheap.top().dist));
         minRadiuses[index] = std::max(minRadiuses[index], nnheap.top().dist);
      }
      nnheap.pop();
   }
   neighborsCount[index] += newNeighborsCount;
#ifdef GENERATE_STATS
   stats.nnCount += newNeighborsCount;
#endif
   if (neighborsCount[index] > n - index || newNeighborsCount == 0)
      shouldFind[index] = false;
   else {
      pq.push(HeapHierarchicalItem(index, SIZE_MAX, minRadiuses[index])); // to be continued...
   }
#ifdef _OPENMP
   omp_unset_lock(&pqwritelock);
#endif
}


void HClustNNbasedSingle::computePrefetch(std::priority_queue< HeapHierarchicalItem > & pq)
{
   // INIT: Pre-fetch a few nearest neighbors for each point
   MESSAGE_2("[%010.3f] prefetching NNs\n", clock()/(float)CLOCKS_PER_SEC);

#ifdef _OPENMP
   omp_set_dynamic(0); /* the runtime will not dynamically adjust the number of threads */
   #pragma omp parallel for schedule(dynamic)
#endif
   for (size_t i=0; i<n; i++)
   {
      if (MASTER_OR_SINGLE_THREAD) Rcpp::checkUserInterrupt(); // may throw an exception, fast op, not thread safe

      getNearestNeighbors(pq, i);

      if (MASTER_OR_SINGLE_THREAD) {
         if (i % 64 == 0) MESSAGE_7("\r             prefetch NN: %d/%d", i, n-1);
      }
   }
   MESSAGE_7("\r             prefetch NN: %d/%d  \n", n-1, n-1);
}


void HClustNNbasedSingle::computeMerge(
      std::priority_queue< HeapHierarchicalItem > & pq,
      HClustResult& res)
{
   MESSAGE_2("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);

   volatile bool go=true;
   volatile size_t i = 0;
#ifdef _OPENMP
   #pragma omp parallel
#endif
   while (go)
   {
#ifdef _OPENMP
      omp_set_lock(&pqwritelock);
#endif
      STOPIFNOT(!pq.empty())
      HeapHierarchicalItem hhi = pq.top();

      if (hhi.index2 == SIZE_MAX) {
         pq.pop();
#ifdef _OPENMP
         omp_unset_lock(&pqwritelock);
#endif
         getNearestNeighbors(pq, hhi.index1);
         continue;
      }

      size_t s1 = ds.find_set(hhi.index1);
      size_t s2 = ds.find_set(hhi.index2);

      if (s1 == s2)
      {
         pq.pop();
#ifdef _OPENMP
         omp_unset_lock(&pqwritelock);
#endif
         continue;
      }

#ifdef _OPENMP
      omp_unset_lock(&pqwritelock); //different threads will be unable to put data into pq without it
      #pragma omp barrier
      #pragma omp single
#endif
      {
         hhi = pq.top(); //it can change, because other threads can push something
         pq.pop();
         s1 = ds.find_set(hhi.index1);
         s2 = ds.find_set(hhi.index2);
         STOPIFNOT(s1 != s2);
         STOPIFNOT(s2 != SIZE_MAX);
         STOPIFNOT(hhi.index1 < hhi.index2);

         res.link(indices[hhi.index1], indices[hhi.index2], hhi.dist);
         ds.link(s1, s2);

         ++i;
         if (i == n-1)
            go = false;/* avoids computing unnecessary nn */
      } // #pragma omp single
      if (MASTER_OR_SINGLE_THREAD) {
         if (i % 512 == 0) MESSAGE_7("\r             merge clusters: %d / %d", i+1, n-1);
         Rcpp::checkUserInterrupt(); // may throw an exception, fast op, not thread safe
      }
   }

   MESSAGE_7("\r             merge clusters: %d / %d  \n", n-1, n-1);
   Rcpp::checkUserInterrupt();
}


HClustResult HClustNNbasedSingle::compute(bool lite)
{
   std::priority_queue< HeapHierarchicalItem > pq;
   // HclustPriorityQueue pq(n);
   HClustResult res(n, distance, lite);

#if VERBOSE >= 5
   distance->getStats().print();
#endif

   prefetch = true;
   computePrefetch(pq);
   prefetch = false;

#if VERBOSE >= 5
   distance->getStats().print();
#endif

   computeMerge(pq, res);

   return res;
}
