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


#include "hclust2_nnbased_nntest_single_approx.h"

using namespace Rcpp;
using namespace std;
using namespace boost;
using namespace grup;


// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustNNbasedNNTestSingleApprox::HClustNNbasedNNTestSingleApprox(Distance* dist, RObject control) :
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
      ds(dist->getObjectCount()),
      prefetch(true),
      vptree(opts, stats, prefetch, ds, n, dist, indices)
{
   NNHeap::setOptions(&opts);

   // starting indices: random permutation of {0,1,...,_n-1}
   for (size_t i=0;i<n;i++)
      indices[i] = i;
   for (size_t i=n-1; i>= 1; i--)
      swap(indices[i], indices[(size_t)(unif_rand()*(i+1))]);
   vptree.build();
}

HClustNNbasedNNTestSingleApprox::~HClustNNbasedNNTestSingleApprox() {

}

void HClustNNbasedNNTestSingleApprox::getNearestNeighbor(size_t index, NumericMatrix mat)
{
   size_t clusterIndex = ds.find_set(index);

#ifdef GENERATE_STATS
#ifdef _OPENMP
#pragma omp atomic
#endif
      ++stats.nnCals;
#endif
   //vptree.zeroNodesVisited();
   NNHeap nnheap;
   vptree.getNearestNeighborsFromMinRadius(index, clusterIndex, minRadiuses[index], nnheap);
   nnheap.fill(nearestNeighbors[index]);

   //Rcout << "Znaleziono " << nearestNeighbors[index].size() << " sasiadow" << endl;
   if (!nearestNeighbors[index].empty())
   {
      int column = 0;
      while (!nearestNeighbors[index].empty()) {
#ifdef GENERATE_STATS
#ifdef _OPENMP
#pragma omp atomic
#endif
         ++stats.nnCount;
#endif
         auto res = nearestNeighbors[index].front();
         nearestNeighbors[index].pop_front();
         //Rcout << "column = " << column << endl;
         //Rcout << "res.dist = " << res.dist << endl;
         if(column < mat.ncol())
            mat(index, column++) = res.dist;
         else
         {
            nearestNeighbors[index].clear();
            break;
         }
      }
   }
}

List HClustNNbasedNNTestSingleApprox::computePrefetch(HclustPriorityQueue& pq)
{
   // INIT: Pre-fetch a few nearest neighbors for each point
   MESSAGE_2("[%010.3f] prefetching NNs\n", clock()/(float)CLOCKS_PER_SEC);

   const int k = opts.minNNPrefetch;
   const vector<int> ograniczeniaLisci = {1, 2, 4, 8, 12, 16, 32, INT_MAX };
   List results(ograniczeniaLisci.size());
   for(size_t i = 0; i< ograniczeniaLisci.size(); ++i)
      results[i] = NumericMatrix(n,k);

#ifdef _OPENMP
   omp_set_dynamic(0); /* the runtime will not dynamically adjust the number of threads */
   omp_lock_t writelock;
   omp_init_lock(&writelock);
   #pragma omp parallel for schedule(dynamic)
#endif
   for(size_t j=0; j < ograniczeniaLisci.size(); ++j)
   {
      opts.nodesVisitedLimit = ograniczeniaLisci[j];
      for (size_t i=0; i<n; i++)
      {
#ifndef _OPENMP
         Rcpp::checkUserInterrupt(); // may throw an exception, fast op, not thread safe
#endif
         NumericMatrix mymat = (NumericMatrix)Rcpp::as<Rcpp::NumericMatrix>(results[j]);
         getNearestNeighbor(i, mymat);
      }
   }
#ifdef _OPENMP
   omp_destroy_lock(&writelock);
#endif
   MESSAGE_7("\r             prefetch NN: %d/%d\n", n-1, n-1);
   return results;
}

List HClustNNbasedNNTestSingleApprox::compute()
{
   HclustPriorityQueue pq(n);
   prefetch = true;
   return computePrefetch(pq);
}
