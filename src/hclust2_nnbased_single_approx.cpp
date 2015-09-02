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


#include "hclust2_nnbased_single_approx.h"

using namespace Rcpp;
using namespace std;
using namespace boost;
using namespace DataStructures;


// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustNNbasedSingleApprox::HClustNNbasedSingleApprox(Distance* dist, RObject control) :
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
      prefetch(true),
      ds(dist->getObjectCount()),
      vptree(opts, stats, prefetch, ds, n, dist, indices)
{
   // starting indices: random permutation of {0,1,...,_n-1}
   for (size_t i=0;i<n;i++)
      indices[i] = i;
   for (size_t i=n-1; i>= 1; i--)
      swap(indices[i], indices[(size_t)(unif_rand()*(i+1))]);
   vptree.build();
}

HClustNNbasedSingleApprox::~HClustNNbasedSingleApprox() {

}

HeapNeighborItem HClustNNbasedSingleApprox::getNearestNeighbor(size_t index, double distMax)
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
      vptree.zeroNodesVisited();
      NNHeap nnheap((prefetch)?opts.maxNNPrefetch:opts.maxNNMerge);
      vptree.getNearestNeighborsFromMinRadius(index, clusterIndex, minRadiuses[index], nnheap);
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

void HClustNNbasedSingleApprox::computePrefetch(HclustPriorityQueue& pq)
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
#ifndef _OPENMP
      Rcpp::checkUserInterrupt(); // may throw an exception, fast op, not thread safe
#endif
      HeapNeighborItem hi=getNearestNeighbor(i);
      if (hi.index != SIZE_MAX)
      {
#if !defined(_OPENMP)
         MESSAGE_7("\r             prefetch NN: %d/%d", i, n-1);
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
   MESSAGE_7("\r             prefetch NN: %d/%d\n", n-1, n-1);
}


void HClustNNbasedSingleApprox::computeMerge(
      HclustPriorityQueue& pq,
      HClustResult& res)
{
   MESSAGE_2("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);

   /*
    * stworz zbior LINKING (czas terazniejszy continous, bo one sie wlasnie teraz linkuja)
    * niech bedzie on pusty
    */

   size_t i = 0;
   while (true)
   {
      //poczatek sekcji krytycznej (sekcja krytyczna powinna byc na dzialania na pq i LINKING)
      HeapHierarchicalItem hhi = pq.top();

      /*
       * tutaj sobie sprawdzaj, czy aby mozna go przetworzyc. jak nie mozna, bo jest w
       * LINKING, to sie zatrzymaj, nie zdejmuj go z kolejki, bo po wykonaniu sie innego linkinga
       * to moze juz nie byc pierwszy element w kolejce. Czekaj na sygnal od innego watku z konca petli
       *
       * tutaj moze miec miejsce taki efekt, ze raz bierzemy element, i mamy s1 i s2, ktore sa przetwarzane,
       * a jak drugi raz wezmiemy ten element, to mozemy dostac s1 i s1, albo np. s1 i s3 czy s3 i s4
       *
       * wersja 1: jak znajdziemy element, ktorego nie mozemy przetworzyc, to stop
       *
       * wersja 2: w takim przypadku wyciagamy kolejny element z kolejki az mozemy go przetworzyc
       *
       * jesli mamy element i wiemy, ze mozemy go przetworzyc, to trzeba zalozyc jakis mutex,
       * zeby zdazyc wrzucic do go LINKING zanim inny watek nam wyjmie kolejny element z kolejki,
       * w ktorym moze byc ten sam cluster
       */

      pq.pop();

      if (hhi.index2 == SIZE_MAX) {
         //wrzuc s1 do LINKING
         //koniec sekcji krytycznej
         HeapNeighborItem hi = getNearestNeighbor(hhi.index1, INFINITY);
         if (isfinite(hi.dist))
            pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));//tutaj wyrzucamy s1 z LINKING (sekcja krytyczna), dajemy sygnal
         continue;
      }

      size_t s1 = ds.find_set(hhi.index1);
      size_t s2 = ds.find_set(hhi.index2);
      /*
       * wrzuc s1 do LINKING
       */
      if (s1 != s2)
      {
         /*
          * wrzuc s2 do LINKING
          */
         //koniec sekcji krytycznej

         Rcpp::checkUserInterrupt(); // may throw an exception, fast op

         res.link(indices[hhi.index1], indices[hhi.index2], hhi.dist);
         ds.link(s1, s2);

         /*
          * wyrzuc s2 z LINKING (?) -- o tyle niepotrzebne, ze juz i tak s2 nie bedzie zwracane przez ds.find_set()
          */

         ++i;
         if (i == n-1) break; /* avoids computing unnecessary nn */
      }
      MESSAGE_7("\r             %d / %d", i+1, n);

      STOPIFNOT(hhi.index1 < hhi.index2);

      HeapNeighborItem hi=getNearestNeighbor(hhi.index1, pq.top().dist);
      //poczatek sekcji krytycznej
      STOPIFNOT(hhi.index1 < hi.index);
      if (isfinite(hi.dist))
         pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));

      /*
       * wyrzuc s1 z LINKING
       * rzuc jakims sygnalem, ktory zwolni watek z poczatku petli, aby wzial element
       * z kolejki i sprawdzil, czy moze go przetworzyc
       */
      //koniec sekcji krytycznej
   }
   MESSAGE_7("\r             %d / %d\n", n, n);
   Rcpp::checkUserInterrupt();
}

HClustResult HClustNNbasedSingleApprox::compute()
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
