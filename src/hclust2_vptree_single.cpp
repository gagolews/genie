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
HClustVpTreeSingle::HClustVpTreeSingle(Distance* dist, RObject control) :
   HClustNNbasedSingle(dist, control),
   _root(NULL),
   distances(vector<double>(_n))
{
   MESSAGE_2("[%010.3f] building vp-tree\n", clock()/(float)CLOCKS_PER_SEC);

   _root = buildFromPoints(0, _n);
}


HClustVpTreeSingle::~HClustVpTreeSingle() {
//   MESSAGE_2("[%010.3f] destroying vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
   if(_root) delete _root;
}


size_t HClustVpTreeSingle::chooseNewVantagePoint(size_t left, size_t right)
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


HClustVpTreeSingleNode* HClustVpTreeSingle::buildFromPoints(size_t left,
   size_t right)
{
#ifdef GENERATE_STATS
   ++stats.nodeCount;
#endif
   if (right - left <= opts.maxLeavesElems)
   {
   #ifdef GENERATE_STATS
      ++stats.leafCount;
   #endif
      HClustVpTreeSingleNode* leaf = new HClustVpTreeSingleNode(left, right);
      // std::sort(_indices.begin()+left, _indices.begin()+right, comparer_gt);
      // leaf->maxindex = _indices[left];
      // leaf->maxindex = _indices[left];
      // for (size_t i=left+1; i<right; ++i)
         // if (_indices[i] > leaf->maxindex)
            // leaf->maxindex = _indices[i];
      // for (size_t i=left; i<right; ++i)
         // _indicesinv[_indices[i]] = i;
      leaf->maxindex = right-1;
      return leaf;
   }

   size_t vpi_idx = chooseNewVantagePoint(left, right);
   std::swap(_indices[left], _indices[vpi_idx]);
   size_t vpi = _indices[left];
   // _indicesinv[vpi] = left;
   size_t median = (right + left) / 2;

   for (size_t i=left+1; i<right; ++i)
      distances[_indices[i]] = (*_distance)(vpi, _indices[i]);

   // std::sort(_indices.begin()+left+1, _indices.begin()+right, DistanceComparatorCached(&distances));
   std::nth_element(_indices.begin()+left+1, _indices.begin() + median, _indices.begin()+right, DistanceComparatorCached(&distances));

// slower -- computes some distances > 1 time
//    std::nth_element(_indices.begin() + left + 1, _indices.begin() + median,  _indices.begin() + right,
//                     DistanceComparator(vpi, _distance));
//    HClustVpTreeSingleNode* node = new HClustVpTreeSingleNode(vpi, left, left+1, (*_distance)(vpi, _indices[median]));

   HClustVpTreeSingleNode* node = new HClustVpTreeSingleNode(vpi, left, left+1, distances[_indices[median]]);

   node->maxindex = left;
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


void HClustVpTreeSingle::getNearestNeighborsFromMinRadiusRecursive(
   HClustVpTreeSingleNode* node, size_t index,
   size_t clusterIndex, double minR, double& maxR, NNHeap& nnheap)
{
   // search within (minR, maxR]
   // if (node == NULL) return; // this should not happen
#ifdef GENERATE_STATS
#ifdef _OPENMP
#pragma omp atomic
#endif
   ++stats.nodeVisit;
#endif

   if (!prefetch && node->sameCluster && clusterIndex == ds.find_set(node->left))
      return;
   if (node->vpindex == SIZE_MAX) { // leaf
      if (!prefetch && !node->sameCluster) {
         size_t commonCluster = ds.find_set(node->left);
         for (size_t i=node->left; i<node->right; ++i) {
            size_t currentCluster = ds.find_set(i);
            if (currentCluster != commonCluster) commonCluster = SIZE_MAX;
            if (currentCluster == clusterIndex) continue;
            if (index >= i) continue;
            double dist2 = (*_distance)(_indices[index], _indices[i]); // the slow part
            if (dist2 > maxR || dist2 <= minR) continue;

            nnheap.insert(i, dist2, maxR);
         }
         if (commonCluster != SIZE_MAX)
            node->sameCluster = true; // set to true (btw, may be true already)
      }
      else /* node->sameCluster */ {
         for (size_t i=node->left; i<node->right; ++i) {
            if (index >= i) continue;
            double dist2 = (*_distance)(_indices[index], _indices[i]); // the slow part
            if (dist2 > maxR || dist2 <= minR) continue;

            nnheap.insert(i, dist2, maxR);
         }
      }
      return; // nothing more to do
   }
   // else // not a leaf

   // first visit the vantage point
   double dist = (*_distance)(_indices[index], _indices[node->left]); // the slow part
   if (index < node->left && dist <= maxR && dist > minR &&
         ds.find_set(node->left) != clusterIndex) {
      nnheap.insert(node->left, dist, maxR);
   }

   if (dist < node->radius) {
      if (dist - maxR <= node->radius && dist + node->radius > minR) {
         if (node->childL && index < node->childL->maxindex)
            getNearestNeighborsFromMinRadiusRecursive(node->childL, index, clusterIndex, minR, maxR, nnheap);
      }

      if (dist + maxR >= node->radius) {
         if (node->childR && index < node->childR->maxindex)
            getNearestNeighborsFromMinRadiusRecursive(node->childR, index, clusterIndex, minR, maxR, nnheap);
      }
   }
   else /* ( dist >= node->radius ) */ {
      if (dist + maxR >= node->radius) {
         if (node->childR && index < node->childR->maxindex)
            getNearestNeighborsFromMinRadiusRecursive(node->childR, index, clusterIndex, minR, maxR, nnheap);
      }

      if (dist - maxR <= node->radius && dist + node->radius > minR) {
         if (node->childL && index < node->childL->maxindex)
            getNearestNeighborsFromMinRadiusRecursive(node->childL, index, clusterIndex, minR, maxR, nnheap);
      }
   }
   if (prefetch || node->sameCluster ||
      (node->childL && !node->childL->sameCluster) ||
      (node->childR && !node->childR->sameCluster)
   ) return;

   // otherwise check if node->sameCluster flag needs updating
   size_t commonCluster = ds.find_set(node->left);
   if (node->childL) {
      size_t currentCluster = ds.find_set(node->childL->left);
      if (currentCluster != commonCluster) return; // not ready yet
   }
   if (node->childR) {
      size_t currentCluster = ds.find_set(node->childR->left);
      if (currentCluster != commonCluster) return; // not ready yet
   }
   node->sameCluster = true;
}


void HClustVpTreeSingle::print(HClustVpTreeSingleNode* n) {
   if (n->childL) {
      Rprintf("\"%llx\" -> \"%llx\" [label=\"L\"];\n",
         (unsigned long long)n, (unsigned long long)(n->childL));
      print(n->childL);
   }
   if (n->childR) {
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


void HClustVpTreeSingle::print() {
   Rprintf("digraph vptree {\n");
   Rprintf("size=\"6,6\";\n");
   Rprintf("node [color=lightblue2, style=filled];");
   print(_root);
   Rprintf("}\n");
}
