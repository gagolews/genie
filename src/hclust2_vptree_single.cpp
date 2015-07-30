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


#include "hclust2_vptree_single.h"



/* Improvement ideas:
 *
 * 1. add custom sort of input objects
 *       useful for Levenshtein distance
 *       long strings should be put at the end
 */


using namespace Rcpp;
using namespace std;
using namespace boost;
using namespace DataStructures;


// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustBiVpTreeSingle::HClustBiVpTreeSingle(Distance* dist, RObject control) :
   opts(control), _root(NULL), _n(dist->getObjectCount()), _distance(dist),
   _indices(dist->getObjectCount()),
   neighborsCount(vector<size_t>(dist->getObjectCount(), 0)),
   minRadiuses(vector<double>(dist->getObjectCount(), -INFINITY)),
   // maxRadiuses(vector<double>(dist->getObjectCount(), INFINITY)),
   shouldFind(vector<bool>(dist->getObjectCount(), true)),
   nearestNeighbors(vector< deque<HeapNeighborItem> >(dist->getObjectCount())),
#ifdef GENERATE_STATS
   stats(HClustBiVpTreeStats()),
#endif
   ds(dist->getObjectCount()),
   heap(std::priority_queue<HeapNeighborItem>())
{
#if VERBOSE > 5
   Rprintf("[%010.3f] building vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
#endif

   // starting _indices: random permutation of {0,1,...,_n-1}
   for (size_t i=0;i<_n;i++)
      _indices[i] = i;
   for (size_t i=_n-1; i>= 1; i--)
      swap(_indices[i], _indices[(size_t)(unif_rand()*(i+1))]);

   _root = buildFromPoints(0, _n);
}


HClustBiVpTreeSingle::~HClustBiVpTreeSingle() {
// #if VERBOSE > 5
//       Rprintf("[%010.3f] destroying vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
// #endif
   if(_root) delete _root;
}


size_t HClustBiVpTreeSingle::chooseNewVantagePoint(size_t left, size_t right)
{
   if (opts.vpSelectScheme == 1) {
      // idea by Yianilos (original vp-tree paper)
      if (left + opts.vpSelectCand + opts.vpSelectTest > right)
         return left;

      // randomize:
      for (size_t i=left; i<left+opts.vpSelectCand+opts.vpSelectTest; ++i)
         std::swap(_indices[i], _indices[i+(size_t)(unif_rand()*(right-i))]);

      // maximize variance
      size_t bestIndex = -1;
      double bestSigma = -INFINITY;
      for (size_t i=left; i<left+opts.vpSelectCand; i++) {
         accumulators::accumulator_set< double,
            accumulators::features<accumulators::tag::variance> > acc;
         for (size_t j = left+opts.vpSelectCand; j < left+opts.vpSelectCand+opts.vpSelectTest; ++j)
            acc( (*_distance)( _indices[i], _indices[j] ) );
         double curSigma = accumulators::variance(acc);
         if (curSigma > bestSigma) {
            bestSigma = curSigma;
            bestIndex = i;
         }
      }

      return bestIndex;
   }
   else if (opts.vpSelectScheme == 2) {
      // idea by T. Bozkaya and M. Ozsoyoglu, "Indexing large metric spaces
      //      for similarity search queries"

      // randomize:
      std::swap(_indices[left], _indices[left+(size_t)(unif_rand()*(right-left))]);

      // which one maximizes dist to _indices[left]?
      size_t bestIndex = left;
      double bestDist  = 0.0;
      for (size_t i=left+1; i<right; ++i) {
         double curDist = (*_distance)(_indices[left], _indices[i]);
         if (curDist > bestDist) {
            bestDist = curDist;
            bestIndex = i;
         }
      }
   //       for (size_t i=left+2; i<right; ++i) {
   //          double curDist = (*_distance)(_indices[left+1], _indices[i]);
   //          if (curDist > bestDist) {
   //             bestDist = curDist;
   //             bestIndex = i;
   //          }
   //       }
   //       for (size_t i=left+3; i<right; ++i) {
   //          double curDist = (*_distance)(_indices[left+2], _indices[i]);
   //          if (curDist > bestDist) {
   //             bestDist = curDist;
   //             bestIndex = i;
   //          }
   //       }
      return bestIndex;
   }
   else {
      // return random index
      // don'use left one (even if sample seems to be randomized already,
      // vp in subtrees is already on the left...)
      return left+(size_t)(unif_rand()*(right-left));
   }
}

bool comparer_gt(int i,int j) { return (i>j); }


HClustBiVpTreeSingleNode* HClustBiVpTreeSingle::buildFromPoints(size_t left, size_t right)
{
#ifdef GENERATE_STATS
   ++stats.nodeCount;
#endif
   if (right - left <= opts.maxLeavesElems)
   {
   #ifdef GENERATE_STATS
      ++stats.leafCount;
   #endif
      HClustBiVpTreeSingleNode* leaf = new HClustBiVpTreeSingleNode(left, right);
      std::sort(_indices.begin()+left, _indices.begin()+right, comparer_gt);
      leaf->maxindex = _indices[left];
      // leaf->maxindex = _indices[left];
      // for (size_t i=left+1; i<right; ++i)
         // if (_indices[i] > leaf->maxindex)
            // leaf->maxindex = _indices[i];
      return leaf;
   }

   size_t vpi_idx = chooseNewVantagePoint(left, right);
   std::swap(_indices[left], _indices[vpi_idx]);
   size_t vpi = _indices[left];

   size_t median = (right + left) / 2;
   std::nth_element(_indices.begin() + left + 1, _indices.begin() + median,  _indices.begin() + right,
                    DistanceComparator(vpi, _distance));

   HClustBiVpTreeSingleNode* node = new HClustBiVpTreeSingleNode(vpi, left, left+1, (*_distance)(vpi, _indices[median]));

   node->maxindex = vpi;
   if (median - left > 0) { // don't include vpi
      node->childL = buildFromPoints(left+1, median+1);
      if (node->childL->maxindex > node->maxindex)
         node->maxindex = node->childL->maxindex;
   }
   if (right - median - 1 > 0) {
      node->childR = buildFromPoints(median+1, right);
      if (node->childR->maxindex > node->maxindex)
         node->maxindex = node->childR->maxindex;
   }

   return node;
}


void HClustBiVpTreeSingle::getNearestNeighborsFromMinRadiusRecursive(
   HClustBiVpTreeSingleNode* node, size_t index,
   size_t clusterIndex, double minR, double& maxR)
{
   // search within (minR, maxR]
   // if (node == NULL) return; // this should not happen
#ifdef GENERATE_STATS
   ++stats.nodeVisit;
#endif

   if (!prefetch && node->sameCluster && clusterIndex == ds.find_set(_indices[node->left]))
      return;

   if (node->vpindex == SIZE_MAX) { // leaf
      if (!prefetch && !node->sameCluster) {
         size_t commonCluster = ds.find_set(_indices[node->left]);
         for (size_t i=node->left; i<node->right; ++i) {
            size_t currentCluster = ds.find_set(_indices[i]);
            if (currentCluster != commonCluster) commonCluster = SIZE_MAX;
            if (currentCluster == clusterIndex) continue;
            if (index >= _indices[i]) continue;
            double dist2 = (*_distance)(index, _indices[i]); // the slow part
            if (dist2 > maxR || dist2 <= minR) continue;

            if (heap.size() >= opts.maxNNPrefetch && dist2 < maxR) {
               while (!heap.empty() && heap.top().dist == maxR) {
                  heap.pop();
               }
            }
            heap.push( HeapNeighborItem(_indices[i], dist2) );
            if (heap.size() >= opts.maxNNPrefetch) maxR = heap.top().dist;
         }
         if (commonCluster != SIZE_MAX)
            node->sameCluster = true; // set to true (btw, may be true already)
      }
      else /* node->sameCluster */ {
         for (size_t i=node->left; i<node->right; ++i) {
            if (index >= _indices[i]) break; // indices are sorted
            double dist2 = (*_distance)(index, _indices[i]); // the slow part
            if (dist2 > maxR || dist2 <= minR) continue;

            if (heap.size() >= opts.maxNNPrefetch && dist2 < maxR) {
               while (!heap.empty() && heap.top().dist == maxR) {
                  heap.pop();
               }
            }
            heap.push( HeapNeighborItem(_indices[i], dist2) );
            if (heap.size() >= opts.maxNNPrefetch) maxR = heap.top().dist;
         }
      }
      return; // nothing more to do
   }
   // else // not a leaf

   // first visit the vantage point
   double dist = (*_distance)(node->vpindex, index); // the slow part
   if (ds.find_set(node->vpindex) != clusterIndex && index < node->vpindex) {
      if (dist <= maxR && dist > minR) {
         if (heap.size() >= opts.maxNNPrefetch && dist < maxR) {
            while (!heap.empty() && heap.top().dist == maxR) {
               heap.pop();
            }
         }
         heap.push( HeapNeighborItem(node->vpindex, dist) );
         if (heap.size() >= opts.maxNNPrefetch) maxR = heap.top().dist;
      }
   }

   if (dist < node->radius) {
      if (dist - maxR <= node->radius && dist + node->radius > minR) {
         if (node->childL && index < node->childL->maxindex)
            getNearestNeighborsFromMinRadiusRecursive(node->childL, index, clusterIndex, minR, maxR);
      }

      if (dist + maxR >= node->radius) {
         if (node->childR && index < node->childR->maxindex)
            getNearestNeighborsFromMinRadiusRecursive(node->childR, index, clusterIndex, minR, maxR);
      }
   }
   else /* ( dist >= node->radius ) */ {
      if (dist + maxR >= node->radius) {
         if (node->childR && index < node->childR->maxindex)
            getNearestNeighborsFromMinRadiusRecursive(node->childR, index, clusterIndex, minR, maxR);
      }

      if (dist - maxR <= node->radius && dist + node->radius > minR) {
         if (node->childL && index < node->childL->maxindex)
            getNearestNeighborsFromMinRadiusRecursive(node->childL, index, clusterIndex, minR, maxR);
      }
   }

   if (prefetch || node->sameCluster ||
      (node->childL && !node->childL->sameCluster) ||
      (node->childR && !node->childR->sameCluster)
   ) return;

   // otherwise check if node->sameCluster flag needs updating
   size_t commonCluster = ds.find_set(node->vpindex);
   if (node->childL) {
      size_t currentCluster = ds.find_set(_indices[node->childL->left]);
      if (currentCluster != commonCluster) return; // not ready yet
   }
   if (node->childR) {
      size_t currentCluster = ds.find_set(_indices[node->childR->left]);
      if (currentCluster != commonCluster) return; // not ready yet
   }
   node->sameCluster = true;
}


HeapNeighborItem HClustBiVpTreeSingle::getNearestNeighbor(size_t index)
{
   if (shouldFind[index] && nearestNeighbors[index].empty())
   {
      size_t clusterIndex = ds.find_set(index);
      double _tau = INFINITY;//maxRadiuses[index];

#ifdef GENERATE_STATS
      ++stats.nnCals;
#endif
      getNearestNeighborsFromMinRadiusRecursive(_root, index, clusterIndex, minRadiuses[index], _tau);
      while (!heap.empty()) {
         nearestNeighbors[index].push_front(heap.top());
         heap.pop();
      }
      // maxRadiuses[index] = INFINITY;
      size_t newNeighborsCount = nearestNeighbors[index].size();

      neighborsCount[index] += newNeighborsCount;
      if (neighborsCount[index] > _n - index || newNeighborsCount == 0)
         shouldFind[index] = false;

      if (newNeighborsCount > 0)
         minRadiuses[index] = nearestNeighbors[index].back().dist;
   }

   if (!nearestNeighbors[index].empty())
   {
#ifdef GENERATE_STATS
      ++stats.nnCount;
#endif
      auto res = nearestNeighbors[index].front();
      nearestNeighbors[index].pop_front();
      return res;
   }
   else
   {
      return HeapNeighborItem(SIZE_MAX,-INFINITY);
   }
}


NumericMatrix HClustBiVpTreeSingle::compute()
{
   NumericMatrix ret(_n-1, 2);
   priority_queue<HeapHierarchicalItem> pq;

   // INIT: Pre-fetch a few nearest neighbors for each point
#if VERBOSE > 5
   Rprintf("[%010.3f] prefetching NNs\n", clock()/(float)CLOCKS_PER_SEC);
#endif

   prefetch = true;
   for (size_t i=0; i<_n; i++)
   {
#if VERBOSE > 7
      if (i % 1024 == 0) Rprintf("\r             prefetch NN: %d/%d", i, _n-1);
#endif
      Rcpp::checkUserInterrupt(); // may throw an exception, fast op
      HeapNeighborItem hi=getNearestNeighbor(i);

      if (hi.index != SIZE_MAX)
      {
         pq.push(HeapHierarchicalItem(i, hi.index, hi.dist));
      }
   }
#if VERBOSE > 7
   Rprintf("\r             prefetch NN: %d/%d\n", _n-1, _n-1);
#endif
#if VERBOSE > 5
   Rprintf("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);
#endif

   prefetch = false;
   size_t i = 0;
   while(i < _n - 1)
   {
      //Rcout << "iteracja " << i << endl;
      //Rcout << "pq size = " << pq.size()<< endl;
      HeapHierarchicalItem hhi = pq.top();
      pq.pop();

      size_t s1 = ds.find_set(hhi.index1);
      size_t s2 = ds.find_set(hhi.index2);
      if (s1 != s2)
      {
         Rcpp::checkUserInterrupt(); // may throw an exception, fast op

         ret(i,0)=(double)hhi.index1;
         ret(i,1)=(double)hhi.index2;
         ++i;
         ds.link(s1, s2);
      }
#if VERBOSE > 7
      if (i % 1024 == 0) Rprintf("\r             %d / %d", i+1, _n);
#endif

      // ASSERT: hhi.index1 < hhi.index2
      HeapNeighborItem hi=getNearestNeighbor(hhi.index1);
      if (hi.index != SIZE_MAX)
         pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));
   }
#if VERBOSE > 7
   Rprintf("\r             %d / %d\n", _n, _n);
#endif
   Rcpp::checkUserInterrupt();

#if VERBOSE > 5
   Rprintf("[%010.3f] generating output matrix\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   MergeMatrixGenerator mmg(ret.nrow());
   return mmg.generateMergeMatrix(ret);
}


void HClustBiVpTreeSingle::print(HClustBiVpTreeSingleNode* n) {
   if (n->childL) {
      Rprintf("\"%llx\" -> \"%llx\" [label=\"L\"];\n",
         (unsigned long long)n, (unsigned long long)(n->childL));
      print(n->childL);
   }
   if (n->childL) {
      Rprintf("\"%llx\" -> \"%llx\" [label=\"R\"];\n",
         (unsigned long long)n, (unsigned long long)(n->childR));
      print(n->childR);
   }

   if (n->vpindex == SIZE_MAX) {
      for (size_t i=n->left; i<n->right; ++i)
         Rprintf("\"%llx\" -> \"%llu\" [arrowhead = diamond];\n", (unsigned long long)n, (unsigned long long)_indices[i]+1);
   }
   else {
      Rprintf("\"%llx\" [label=\"(%llu, %g)\"];\n", (unsigned long long)n, (unsigned long long)n->vpindex+1, n->radius);
   }
}


void HClustBiVpTreeSingle::print() {
   Rprintf("digraph vptree {\n");
   Rprintf("size=\"6,6\";\n");
   Rprintf("node [color=lightblue2, style=filled];");
   print(_root);
   Rprintf("}\n");
}


// [[Rcpp::export(".hclust2_single")]]
RObject hclust2_single(RObject distance, RObject objects, RObject control=R_NilValue) {
#if VERBOSE > 5
   Rprintf("[%010.3f] starting timer\n", clock()/(double)CLOCKS_PER_SEC);
#endif
   RObject result(R_NilValue);
   DataStructures::Distance* dist = DataStructures::Distance::createDistance(distance, objects);

   try {
      /* Rcpp::checkUserInterrupt(); may throw an exception */
      DataStructures::HClustBiVpTreeSingle hclust(dist, control);
      RObject merge = hclust.compute();
      result = Rcpp::as<RObject>(List::create(
         _["merge"]  = merge,
         _["height"] = R_NilValue,
         _["order"]  = R_NilValue,
         _["labels"] = R_NilValue,
         _["call"]   = R_NilValue,
         _["method"] = "single",
         _["dist.method"] = R_NilValue,
         _["stats"] = List::create(
            _["vptree"] = hclust.getStats().toR(),
            _["distance"] = dist->getStats().toR()
         ),
         _["control"] = List::create(
            _["vptree"] = hclust.getOptions().toR()
         )
      ));
      result.attr("class") = "hclust";
      //hclust.print();
   }
   catch(...) {

   }

   if (dist) delete dist;
#if VERBOSE > 5
   Rprintf("[%010.3f] done\n", clock()/(double)CLOCKS_PER_SEC);
#endif
   if (Rf_isNull(result)) stop("stopping on error or explicit user interrupt");
   return result;
}
