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
      n(dist->getObjectCount()),
      distance(dist),
      indices(dist->getObjectCount()),
      neighborsCount(dist->getObjectCount(), 0),
      minRadiuses(dist->getObjectCount(), -INFINITY),
      shouldFind(dist->getObjectCount(), true),
      nearestNeighbors(dist->getObjectCount()),
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
      if (neighborsCount[index] > n - index || newNeighborsCount == 0)
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


void HClustNNbasedSingle::computePrefetch(HclustPriorityQueue& pq)
{
   // INIT: Pre-fetch a few nearest neighbors for each point
   MESSAGE_2("[%010.3f] prefetching NNs\n", clock()/(float)CLOCKS_PER_SEC);

#ifdef _OPENMP
   omp_set_dynamic(0); /* the runtime will not dynamically adjust the number of threads */
   omp_lock_t writelock;
   omp_init_lock(&writelock);
   #pragma omp parallel for schedule(dynamic)
#endif
   for (size_t i=0; i<n; i++)
   {
      if (MASTER_OR_SINGLE_THREAD) Rcpp::checkUserInterrupt(); // may throw an exception, fast op, not thread safe

      HeapNeighborItem hi=getNearestNeighbor(i);
      if (hi.index != SIZE_MAX)
      {
         if (MASTER_OR_SINGLE_THREAD) {
            if (i % 64 == 0) MESSAGE_7("\r             prefetch NN: %d/%d", i, n-1);
         }
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
   MESSAGE_7("\r             prefetch NN: %d/%d  \n", n-1, n-1);
}


void HClustNNbasedSingle::computeMerge(
      HclustPriorityQueue& pq,
      HClustResult& res)
{
   MESSAGE_2("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);

#ifdef _OPENMP
   int threadMerge;
   omp_set_dynamic(0); /* the runtime will not dynamically adjust the number of threads */
   omp_lock_t writelock; //critical section for pq
   omp_init_lock(&writelock);
#endif

   volatile bool go=true;
   volatile size_t i = 0;
#ifdef _OPENMP
   #pragma omp parallel shared(go, i, pq, res, threadMerge)
#endif
   while (go)
   {
#ifdef _OPENMP
      omp_set_lock(&writelock);
#endif
      HeapHierarchicalItem hhi = pq.top();

      if (hhi.index2 == SIZE_MAX) {
         pq.pop();
#ifdef _OPENMP
         omp_unset_lock(&writelock);
#endif
         HeapNeighborItem hi = getNearestNeighbor(hhi.index1, INFINITY);
         if (isfinite(hi.dist))
         {
#ifdef _OPENMP
            omp_set_lock(&writelock);
#endif
            pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));
#ifdef _OPENMP
            omp_unset_lock(&writelock);
#endif
         }
         continue;
      }

      size_t s1 = ds.find_set(hhi.index1);
      size_t s2 = ds.find_set(hhi.index2);

      if (s1==s2)
      {
         pq.pop();
#ifdef _OPENMP
         omp_unset_lock(&writelock);
#endif
         STOPIFNOT(hhi.index1 < hhi.index2);
         HeapNeighborItem hi=getNearestNeighbor(hhi.index1, pq.top().dist);
         STOPIFNOT(hhi.index1 < hi.index);
         if (isfinite(hi.dist))
         {
#ifdef _OPENMP
            omp_set_lock(&writelock);
#endif
            pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));
#ifdef _OPENMP
            omp_unset_lock(&writelock);
#endif
         }
         continue;
      }

#ifdef _OPENMP
      omp_unset_lock(&writelock); //different threads will be unable to put data into pq without it
      #pragma omp barrier
      #pragma omp single
#endif
      {
#ifdef _OPENMP
         threadMerge = omp_get_thread_num();
#endif
         hhi = pq.top(); //it can change, because other threads can push something
         pq.pop();
         s1 = ds.find_set(hhi.index1);
         s2 = ds.find_set(hhi.index2);
         STOPIFNOT(s1 != s2);
         STOPIFNOT(hhi.index1 < hhi.index2);

         res.link(indices[hhi.index1], indices[hhi.index2], hhi.dist);
         ds.link(s1, s2);

         ++i;
         if (i == n-1)
            go = false;/* avoids computing unnecessary nn */

      }
#ifdef _OPENMP

      if(threadMerge == omp_get_thread_num())
#endif
      {
         if(go) {
            STOPIFNOT(hhi.index1 < hhi.index2);
            HeapNeighborItem hi=getNearestNeighbor(hhi.index1, pq.top().dist);
            STOPIFNOT(hhi.index1 < hi.index);
            if (isfinite(hi.dist))
            {
#ifdef _OPENMP
               omp_set_lock(&writelock);
#endif
               pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));
#ifdef _OPENMP
               omp_unset_lock(&writelock);
#endif
            }
         }
      } // #pragma omp single
      if (MASTER_OR_SINGLE_THREAD) {
         if (i % 512 == 0) MESSAGE_7("\r             merge clusters: %d / %d", i+1, n-1);
         Rcpp::checkUserInterrupt(); // may throw an exception, fast op, not thread safe
      }
   }
#ifdef _OPENMP
   omp_destroy_lock(&writelock);
#endif
   MESSAGE_7("\r             merge clusters: %d / %d  \n", n-1, n-1);
   Rcpp::checkUserInterrupt();
}


HClustResult HClustNNbasedSingle::compute()
{
   HclustPriorityQueue pq(n);
   HClustResult res(n, distance);

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
