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


#include "hclust2_nnbased_gini.h"

using namespace Rcpp;
using namespace std;
using namespace boost;
using namespace grup;


// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustNNbasedGini::HClustNNbasedGini(Distance* dist, RObject control) :
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
   if (opts.thresholdGini < 1.0)
      symmetric = true;
   else
      symmetric = false;

   // starting indices: random permutation of {0,1,...,_n-1}
   for (size_t i=0;i<n;i++)
      indices[i] = i;
   for (size_t i=n-1; i>= 1; i--)
      swap(indices[i], indices[(size_t)(unif_rand()*(i+1))]);
}


HClustNNbasedGini::~HClustNNbasedGini() {

}



HeapNeighborItem HClustNNbasedGini::getNearestNeighbor(size_t index, double distMax)
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


void HClustNNbasedGini::prefetchNNsSymmetric()
{
   std::vector<NNHeap> nnheaps(n);
   std::vector<double> maxR(n, INFINITY);

#ifdef _OPENMP
   std::vector<omp_lock_t> writelocks;
   for (size_t i=0; i<n; ++i) {
      omp_init_lock(&writelocks[i]);
   }
   #pragma omp parallel for schedule(dynamic)
#endif
   for (size_t i=0; i<n; ++i) {
      for (size_t j=i+1; j<n; ++j) {
         double dist2 = (*distance)(indices[i], indices[j]); // the slow part

         OPENMP_ONLY(omp_set_lock(&writelocks[i]))
         nnheaps[i].insert(j, dist2, maxR[i]);
         OPENMP_ONLY(omp_unset_lock(&writelocks[i]))

         OPENMP_ONLY(omp_set_lock(&writelocks[j]))
         nnheaps[j].insert(i, dist2, maxR[j]);
         OPENMP_ONLY(omp_unset_lock(&writelocks[j]))
      }
   }

#ifdef _OPENMP
   #pragma omp parallel for schedule(dynamic)
#endif
   for (size_t i=0; i<n; ++i) {
      nnheaps[i].fill(nearestNeighbors[i]);
   }

#ifdef _OPENMP
   for (size_t i=0; i<n; ++i) {
      omp_destroy_lock(&writelocks[i]);
   }
#endif

}


void HClustNNbasedGini::computePrefetch(HclustPriorityQueue& pq)
{
   // INIT: Pre-fetch a few nearest neighbors for each point
   MESSAGE_2("[%010.3f] prefetching NNs\n", clock()/(float)CLOCKS_PER_SEC);

   if (symmetric && !opts.useVpTree)
      prefetchNNsSymmetric();
   else
      ; /* use vp-tree or something else */

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
         OPENMP_ONLY(omp_set_lock(&writelock))
         pq.push(HeapHierarchicalItem(i, hi.index, hi.dist));
         OPENMP_ONLY(omp_unset_lock(&writelock))
      }
   }
#ifdef _OPENMP
   omp_destroy_lock(&writelock);
#endif
   MESSAGE_7("\r             prefetch NN: %d/%d  \n", n-1, n-1);
}



// double computeGini(PhatDisjointSets& ds, size_t n) {
//    // compute the Gini  coefficient
//    std::size_t start = ds.find_set(0);
//    double gini = 0.0;
//    std::size_t curi = start;
//    do {
//       std::size_t curj = ds.getClusterNext(curi);
//       while (curj != start) {
//          gini += std::fabs((double)ds.getClusterSize(curi)-ds.getClusterSize(curj));
//          curj = ds.getClusterNext(curj);
//       }
//       curi = ds.getClusterNext(curi);
//    } while (curi != start);
//    gini /= (n-1)*(double)ds.getClusterCount();
//    return gini;
// }


void HClustNNbasedGini::linkAndRecomputeGini(double& lastGini, size_t s1, size_t s2)
{
   double size1 = ds.getClusterSize(s1);
   double size2 = ds.getClusterSize(s2);
   lastGini *= (n)*(double)(ds.getClusterCount()-1);
   std::size_t curi = s1;
   do {
      double curs = ds.getClusterSize(curi);
      lastGini -= std::fabs(curs-size1);
      lastGini -= std::fabs(curs-size2);
      lastGini += std::fabs(curs-size1-size2);
      curi = ds.getClusterNext(curi);
   } while (curi != s1);
   lastGini += std::fabs(size2-size1);
   lastGini -= std::fabs(size2-size1-size2);
   lastGini -= std::fabs(size1-size1-size2);

   ds.link(s1, s2);

   lastGini /= (n)*(double)(ds.getClusterCount()-1);
}


// std::size_t getMinClusterSize(PhatDisjointSets& ds) {
//    // static double ttot = 0.0;
//    // double t0 = clock()/(float)CLOCKS_PER_SEC;
//    std::size_t start = ds.find_set(0);
//    std::size_t minsize = ds.getClusterSize(start);
//    size_t curi = ds.getClusterNext(start);
//    while (curi != start) {
//       minsize = min(minsize, ds.getClusterSize(curi));
//       curi = ds.getClusterNext(curi);
//    }
//    // ttot += (clock()/(float)CLOCKS_PER_SEC)-t0;
//    // if (ds.getClusterCount() == 2) cerr << endl << ttot << endl;
//    return minsize;
// }




void HClustNNbasedGini::computeMerge(
      HclustPriorityQueue& pq,
      HClustResult& res)
{
   MESSAGE_2("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);

   #ifdef _OPENMP
   omp_set_dynamic(0); /* the runtime will not dynamically adjust the number of threads */
   omp_lock_t writelock; //critical section for pq
   omp_init_lock(&writelock);
   #endif

   double lastGini = 0.0;
   bool go = true;
   size_t i = 0;
   std::size_t minsize = 1;
   std::deque<HeapHierarchicalItem> pq_cache;
   #ifdef _OPENMP
   #pragma omp parallel shared(go, i, pq, res, lastGini, pq_cache, minsize)
   #endif
   while (go)
   {
      STOPIFNOT(lastGini >= 0 && lastGini <= 1)
      OPENMP_ONLY(omp_set_lock(&writelock))
      // pq may be empty if we have all the elements in pq_cache
      if (!pq.empty()) {
         HeapHierarchicalItem hhi = pq.top();
         STOPIFNOT(isfinite(hhi.dist));
         size_t s1 = ds.find_set(hhi.index1);
         size_t s2 = (hhi.index2 == SIZE_MAX)?s1:ds.find_set(hhi.index2);
         if (s1 == s2) {
            // Two cases are possible here:
            //   a) (hhi.index2 == SIZE_MAX) <-- this was a "fake" PQ element
            //   (used an estimate for the lower bound of the distance to the
            //   next nearest neighbor);
            //   b) we just got two elements in the same cluster;
            // So now it's time to fetch the next "real" NN of hhi.index1
            // the one we get will surely be s.t. s1 != s2
            // the writelock is still in ON
            pq.pop();
            OPENMP_ONLY(omp_unset_lock(&writelock))
            HeapNeighborItem hi = getNearestNeighbor(hhi.index1, INFINITY);
            if (isfinite(hi.dist)) {
               OPENMP_ONLY(omp_set_lock(&writelock))
               pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));
               OPENMP_ONLY(omp_unset_lock(&writelock))
            }
            continue;
         }

         STOPIFNOT(s1 != s2)
         // if lastGini is above thresholdGini, we are only interested in
         // pq elems that are from clusters of size equal to minsize
         if (lastGini > opts.thresholdGini &&
               ds.getClusterSize(s1) > minsize &&
               ds.getClusterSize(s2) > minsize) {
            // the writelock is still in ON
            pq_cache.push_back(hhi);
            pq.pop();
            OPENMP_ONLY(omp_unset_lock(&writelock))
            continue;
         }
      }
      // THREAD BARRIER FOLLOWS
      // all the threads no longer work
      // if we got here, then either we got a merge candidate,
      // or pq is empty

      #ifdef _OPENMP
      omp_unset_lock(&writelock);
      #pragma omp barrier
      #pragma omp single
      #endif
      {
         std::size_t lastminsize = minsize;
         if (!pq.empty()) {
            HeapHierarchicalItem hhi = pq.top();
            pq.pop();
            STOPIFNOT(hhi.index2 != SIZE_MAX)
            size_t s1 = ds.find_set(hhi.index1);
            size_t s2 = ds.find_set(hhi.index2);
            STOPIFNOT(s1 != s2)
            STOPIFNOT(lastGini <= opts.thresholdGini ||
               (ds.getClusterSize(s1) == minsize || ds.getClusterSize(s2) == minsize))

            res.link(indices[hhi.index1], indices[hhi.index2],
               (lastGini <= opts.thresholdGini)?hhi.dist:-hhi.dist);
            linkAndRecomputeGini(lastGini, s1, s2);
            minsize = ds.getMinClusterSize();

            if (++i == n-1)
               go = false;
            pq.push(HeapHierarchicalItem(hhi.index1, SIZE_MAX, hhi.dist));
         }

         if (go && (pq.empty() || lastGini <= opts.thresholdGini || minsize != lastminsize)) {
            if (pq_cache.size() > 5) pq.reset(); // will call make_heap on next top()
            while (!pq_cache.empty()) {
               pq.push(pq_cache.back());
               pq_cache.pop_back();
            }
         }
      } // END OMP BARRIER SINGLE

      if (MASTER_OR_SINGLE_THREAD) {
         if (i % 512 == 0) MESSAGE_7("\r             merge clusters: %d / %d", i+1, n-1);
         Rcpp::checkUserInterrupt(); // may throw an exception, fast op, not thread safe
      }
   } // END WHILE

   #ifdef _OPENMP
   omp_destroy_lock(&writelock);
   #endif
   MESSAGE_7("\r             merge clusters: %d / %d  \n", n-1, n-1);
   Rcpp::checkUserInterrupt();
}


HClustResult HClustNNbasedGini::compute()
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
