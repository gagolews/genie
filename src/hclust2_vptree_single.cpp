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


#include "hclust2_vptree_single.h"

using namespace grup;


// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustVpTreeSingle::HClustVpTreeSingle(Distance* dist, HClustOptions* opts) :
      HClustNNbasedSingle(dist, opts),
      root(NULL)
//    visitAll(false)
{
   MESSAGE_2("[%010.3f] building vp-tree\n", clock()/(float)CLOCKS_PER_SEC);

   std::vector<double> distances(n);
   root = buildFromPoints(0, n, distances);
}


HClustVpTreeSingle::~HClustVpTreeSingle() {
//   MESSAGE_2("[%010.3f] destroying vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
   if(root) delete root;
}


size_t HClustVpTreeSingle::chooseNewVantagePoint(size_t left, size_t right)
{
   // if (opts->vpSelectScheme == 1) {
   //    // idea by Yianilos (original vp-tree paper)
   //    if (left + opts->vpSelectCand + opts->vpSelectTest > right)
   //       return left;
   //
   //    // randomize:
   //    for (size_t i=left; i<left+opts->vpSelectCand+opts->vpSelectTest; ++i)
   //       std::swap(indices[i], indices[i+(size_t)(unif_rand()*(right-i))]);
   //
   //    // maximize variance
   //    size_t bestIndex = -1;
   //    double bestSigma = -INFINITY;
   //    for (size_t i=left; i<left+opts->vpSelectCand; i++) {
   //       accumulators::accumulator_set< double,
   //          accumulators::features<accumulators::tag::variance> > acc;
   //       for (size_t j = left+opts->vpSelectCand; j < left+opts->vpSelectCand+opts->vpSelectTest; ++j)
   //          acc( (*distance)( indices[i], indices[j] ) );
   //       double curSigma = accumulators::variance(acc);
   //       if (curSigma > bestSigma) {
   //          bestSigma = curSigma;
   //          bestIndex = i;
   //       }
   //    }
   //
   //    return bestIndex;
   // }
   // else
   if (opts->vpSelectScheme == 2) {
      // idea by T. Bozkaya and M. Ozsoyoglu, "Indexing large metric spaces
      //      for similarity search queries"

      // randomize:
      std::swap(indices[left], indices[left+(size_t)(unif_rand()*(right-left))]);

      // which one maximizes dist to indices[left]?
      size_t bestIndex = left;
      double bestDist  = 0.0;
      for (size_t i=left+1; i<right; ++i) {
         double curDist = (*distance)(indices[left], indices[i]);
         if (curDist > bestDist) {
            bestDist = curDist;
            bestIndex = i;
         }
      }
   //       for (size_t i=left+2; i<right; ++i) {
   //          double curDist = (*distance)(indices[left+1], indices[i]);
   //          if (curDist > bestDist) {
   //             bestDist = curDist;
   //             bestIndex = i;
   //          }
   //       }
   //       for (size_t i=left+3; i<right; ++i) {
   //          double curDist = (*distance)(indices[left+2], indices[i]);
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
   size_t right, std::vector<double>& distances)
{
#ifdef GENERATE_STATS
   ++stats.nodeCount;
#endif
   if (right - left <= opts->maxLeavesElems)
   {
   #ifdef GENERATE_STATS
      ++stats.leafCount;
   #endif
      HClustVpTreeSingleNode* leaf = new HClustVpTreeSingleNode(left, right);
      leaf->maxindex = right-1; // left < right-1
      return leaf;
   }

   size_t vpi_idx = chooseNewVantagePoint(left, right);
   std::swap(indices[left], indices[vpi_idx]);
   size_t vpi = indices[left];
   size_t median = (right + left) / 2;

   for (size_t i=left+1; i<right; ++i)
      distances[indices[i]] = (*distance)(vpi, indices[i]);

   // std::sort(indices.begin()+left+1, indices.begin()+right, DistanceComparatorCached(&distances));
   std::nth_element(indices.begin()+left+1, indices.begin() + median, indices.begin()+right, DistanceComparatorCached(&distances));

// slower -- computes some distances > 1 time
//    std::nth_element(indices.begin() + left + 1, indices.begin() + median,  indices.begin() + right,
//                     DistanceComparator(vpi, distance));
//    HClustVpTreeSingleNode* node = new HClustVpTreeSingleNode(vpi, left, left+1, (*distance)(vpi, indices[median]));

   HClustVpTreeSingleNode* node = new HClustVpTreeSingleNode(vpi, left, left+1, distances[indices[median]]);

   node->maxindex = left;
   if (median - left > 0) { // don't include vpi
      node->childL = buildFromPoints(left+1, median+1, distances);
      if (node->childL->maxindex > node->maxindex)
         node->maxindex = node->childL->maxindex;
   }
   if (right - median - 1 > 0) {
      node->childR = buildFromPoints(median+1, right, distances);
      if (node->childR->maxindex > node->maxindex)
         node->maxindex = node->childR->maxindex;
   }

   return node;
}


void HClustVpTreeSingle::getNearestNeighborsFromMinRadiusRecursiveLeaf(
   HClustVpTreeSingleNode* node, size_t index,
   size_t clusterIndex, double minR, std::priority_queue<double>& bestR, double& maxR, NNHeap& nnheap)
{
   STOPIFNOT(node->vpindex == SIZE_MAX);
   if (!prefetch && !node->sameCluster) {
      size_t commonCluster = ds.find_set(node->left);
      for (size_t i=node->left; i<node->right; ++i) {
         size_t currentCluster = ds.find_set(i);
         if (currentCluster != commonCluster) commonCluster = SIZE_MAX;
         if (currentCluster == clusterIndex) continue;
         if (index >= i) continue;
         double dist2 = (*distance)(indices[index], indices[i]); // the slow part
         if (dist2 > maxR || dist2 <= minR) continue;
         if (dist2 < bestR.top()) { bestR.pop(); bestR.push(dist2); }

         nnheap.insert(i, dist2, maxR);
      }
      if (commonCluster != SIZE_MAX)
         node->sameCluster = true; // set to true (btw, may be true already)
   }
   else /* node->sameCluster */ {
      for (size_t i=node->left; i<node->right; ++i) {
         if (index >= i) continue;
         double dist2 = (*distance)(indices[index], indices[i]); // the slow part
         if (dist2 > maxR || dist2 <= minR) continue;
         if (dist2 < bestR.top()) { bestR.pop(); bestR.push(dist2); }

         nnheap.insert(i, dist2, maxR);
      }
   }
}


void HClustVpTreeSingle::getNearestNeighborsFromMinRadiusRecursiveNonLeaf(
   HClustVpTreeSingleNode* node, size_t index,
   size_t clusterIndex, double minR, std::priority_queue<double>& bestR, double& maxR, NNHeap& nnheap)
{
   STOPIFNOT(node->vpindex != SIZE_MAX);

   // first visit the vantage point
   double dist = (*distance)(indices[index], indices[node->left]); // the slow part
   if (index < node->left && dist <= maxR && dist > minR &&
         ds.find_set(node->left) != clusterIndex) {
      if (dist < bestR.top()) { bestR.pop(); bestR.push(dist); }
      nnheap.insert(node->left, dist, maxR);
   }

//    if (visitAll) {
//       if (node->childL && index < node->childL->maxindex)
//          getNearestNeighborsFromMinRadiusRecursive(node->childL, index, clusterIndex, minR, maxR, nnheap);
//       if (node->childR && index < node->childR->maxindex)
//          getNearestNeighborsFromMinRadiusRecursive(node->childR, index, clusterIndex, minR, maxR, nnheap);
//    }
//    else {
      if (dist < node->radius) {
         if (node->childL && index < node->childL->maxindex && dist + node->radius > minR) {
            // double cutR = dist - node->radius;
            //STOPIFNOT(maxR >= cutR);
            //STOPIFNOT(!(bestR.top() < cutR));
            getNearestNeighborsFromMinRadiusRecursive(node->childL, index, clusterIndex, minR, bestR, maxR, nnheap);
         }

         if (node->childR && index < node->childR->maxindex) {
            double cutR = node->radius - dist;
            if (maxR >= cutR) {
               if (bestR.top() < cutR) {
                  while (!nnheap.empty() && nnheap.top().dist > cutR) {
                     nnheap.pop();
                  }
                  maxR = cutR;
               }
               else
                  getNearestNeighborsFromMinRadiusRecursive(node->childR, index, clusterIndex, minR, bestR, maxR, nnheap);
            }
         }
      }
      else /* ( dist >= node->radius ) */ {
         if (node->childR && index < node->childR->maxindex) {
            // double cutR = node->radius - dist;
            //STOPIFNOT(maxR >= cutR);
            //STOPIFNOT(!(bestR.top() < cutR));
            getNearestNeighborsFromMinRadiusRecursive(node->childR, index, clusterIndex, minR, bestR, maxR, nnheap);
         }

         if (node->childL && index < node->childL->maxindex && dist + node->radius > minR) {
            double cutR = dist - node->radius;
            if (maxR >= cutR) {
               if (bestR.top() < cutR) {
                  while (!nnheap.empty() && nnheap.top().dist > cutR) {
                     nnheap.pop();
                  }
                  maxR = cutR;
               }
               else
                  getNearestNeighborsFromMinRadiusRecursive(node->childL, index, clusterIndex, minR, bestR, maxR, nnheap);
            }
         }
      }
//   }

   updateSameClusterFlag(node);
}


void HClustVpTreeSingle::updateSameClusterFlag(HClustVpTreeSingleNode* node)
{
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


void HClustVpTreeSingle::print(HClustVpTreeSingleNode* node) {
   if (node->childL) {
      Rprintf("\"%llx\" -> \"%llx\" [label=\"L\"];\n",
         (unsigned long long)node, (unsigned long long)(node->childL));
      print(node->childL);
   }
   if (node->childR) {
      Rprintf("\"%llx\" -> \"%llx\" [label=\"R\"];\n",
         (unsigned long long)node, (unsigned long long)(node->childR));
      print(node->childR);
   }

   if (node->vpindex == SIZE_MAX) {
      for (size_t i=node->left; i<node->right; ++i)
         Rprintf("\"%llx\" -> \"%llu\" [arrowhead = diamond];\n", (unsigned long long)node, (unsigned long long)indices[i]+1);
   }
   else {
      Rprintf("\"%llx\" [label=\"(%llu, %g)\"];\n", (unsigned long long)node, (unsigned long long)node->vpindex+1, node->radius);
   }
}


void HClustVpTreeSingle::print() {
   Rprintf("digraph vptree {\n");
   Rprintf("size=\"6,6\";\n");
   Rprintf("node [color=lightblue2, style=filled];");
   print(root);
   Rprintf("}\n");
}
