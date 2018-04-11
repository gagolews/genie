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


#include "hclust2_mstbased_gini.h"
#include "hclust2_vptree_single.h"

using namespace grup;


// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustMSTbasedGini::HClustMSTbasedGini(Distance* dist, HClustOptions* opts) :
      opts(opts),
      n(dist->getObjectCount()),
#ifdef GENERATE_STATS
      stats(),
#endif
      distance(dist)
{

}


HClustMSTbasedGini::~HClustMSTbasedGini()
{
   /* pass */
}


HclustPriorityQueue HClustMSTbasedGini::getMST()
{
   MESSAGE_2("[%010.3f] determing the MST\n", clock()/(float)CLOCKS_PER_SEC);

   HclustPriorityQueue out(n);

   vector<size_t> todo(n-1); // elements which are still not in the spanning tree
   for (size_t k=0; k<n-1; ++k) todo[k] = k+1;
   std::vector<double> Adist(n, INFINITY);
   std::vector<size_t> Afrom(n, SIZE_MAX);

   size_t lastj = 0; // a randomly chosen element :)
   for (size_t i=0; i<n-1; ++i) { // there are n-1 edges in a spanning tree
      size_t bestj = 0;   // always Adist[0] == INFINITY
      size_t bestjpos = 0;

      STOPIFNOT(todo.size() == n-i-1)
      #ifdef _OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (size_t k=0; k<n-i-1; ++k) {
         // the thread-safe part:
         size_t j = todo[k];
         double curdist = (*distance)(lastj,j); // this takes some time...
         if (curdist < Adist[j]) {
            Adist[j] = curdist;
            Afrom[j] = lastj;
         }

         // the thread-unsafe part:
         #ifndef _OPENMP
         if (Adist[bestj] > Adist[j]) {
            bestj = j;
            bestjpos = k;
         }
         #endif
      }

      #ifdef _OPENMP
      // to avoid establishing a (slow!) critical section in
      // the above loop, this fast part is done single-threadedly
      for (size_t k=0; k<n-i-1; ++k) {
         if (Adist[bestj] > Adist[todo[k]]) {
            bestj = todo[k];
            bestjpos = k;
         }
      }
      #endif

      out.push(HeapHierarchicalItem(Afrom[bestj], bestj, Adist[bestj]));
      todo.erase(todo.begin()+bestjpos); // the algorithm is O(n^2) anyway + we need to iterate thru todo sequentially
      lastj = bestj;

      if (i % 512 == 0) MESSAGE_7("\r                    get MST: %d / %d", i, n-1);
      Rcpp::checkUserInterrupt(); // may throw an exception, fast op, not thread safe
   }
   MESSAGE_7("\r                    get MST: %d / %d                   \n", n-1, n-1);

   return out;
}


void HClustMSTbasedGini::linkAndRecomputeGini(PhatDisjointSets& ds, double& lastGini, size_t s1, size_t s2)
{
   // if opts.thresholdGini == 1.0, there's no need to compute the Gini index
   if (opts->thresholdGini < 1.0) {
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
   }

   s1 = ds.link(s1, s2);

   if (opts->thresholdGini < 1.0) {
      lastGini /= (n)*(double)(ds.getClusterCount()-1);
      lastGini = std::min(1.0, std::max(0.0, lastGini)); // avoid numeric inaccuracies
   }
}



HClustResult HClustMSTbasedGini::compute()
{
   HclustPriorityQueue pq;

   if (opts->useVpTree) {
      HClustVpTreeSingle hclust(distance, opts);
      HClustResult res = hclust.compute(/*merge,order not needed*/opts->thresholdGini < 1.0);
      if (opts->thresholdGini >= 1.0) return res;

      Rcpp::NumericMatrix links = res.getLinks();
      Rcpp::NumericVector dist  = res.getHeight();
      STOPIFNOT((size_t)dist.size() == n-1);
      STOPIFNOT((size_t)links.nrow() == n-1);
      pq = HclustPriorityQueue(n);
      for (size_t i=0; i<n-1; ++i) {
         pq.push(HeapHierarchicalItem((size_t)links(i,0), (size_t)links(i,1), (double)dist[i]));
      }
   }
   else {
      pq = getMST();
   }

   HClustResult res(n, distance);
   PhatDisjointSets ds(n);

   MESSAGE_2("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);

   double lastGini = 0.0;
   size_t i = 0;
   std::size_t minsize = 1;
   std::deque<HeapHierarchicalItem> pq_cache;
   while (true)
   {
      STOPIFNOT(!pq.empty())
      HeapHierarchicalItem hhi = pq.top();
      pq.pop();
      size_t s1 = ds.find_set(hhi.index1);
      size_t s2 = ds.find_set(hhi.index2);
      STOPIFNOT(s1 != s2);
      if (lastGini > opts->thresholdGini &&
            ds.getClusterSize(s1) > minsize &&
            ds.getClusterSize(s2) > minsize) {
            // the writelock is still in ON
         pq_cache.push_back(hhi);
         continue;
      }

      std::size_t lastminsize = minsize;
      res.link(hhi.index1, hhi.index2,
         (lastGini <= opts->thresholdGini)?hhi.dist:-hhi.dist);
      linkAndRecomputeGini(ds, lastGini, s1, s2);

      if (opts->thresholdGini < 1.0)
         minsize = ds.getMinClusterSize();

      if (++i == n-1) break;

      if (pq_cache.size() > 0 && (pq.empty() || lastGini <= opts->thresholdGini || minsize != lastminsize))
      {
         if (pq_cache.size() > 5) pq.reset(); // will call make_heap on next top()
         while (!pq_cache.empty()) {
            pq.push(pq_cache.back());
            pq_cache.pop_back();
         }
      }

      if (i % 512 == 0) MESSAGE_7("\r             merge clusters: %d / %d", i+1, n-1);
      Rcpp::checkUserInterrupt(); // may throw an exception, fast op, not thread safe
   } // END WHILE

   MESSAGE_7("\r             merge clusters: %d / %d            \n", n-1, n-1);
   Rcpp::checkUserInterrupt();

   return res;
}
