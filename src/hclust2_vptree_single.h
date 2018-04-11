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

#ifndef __HCLUST2_VPTREE_SINGLE_H
#define __HCLUST2_VPTREE_SINGLE_H



// ************************************************************************

#include "hclust2_nnbased_single.h"


namespace grup
{

struct HClustVpTreeSingleNode
{
   size_t vpindex;
   size_t left;
   size_t right;
   double radius;
   bool sameCluster;
   size_t maxindex;
   HClustVpTreeSingleNode* childL;
   HClustVpTreeSingleNode* childR;

   HClustVpTreeSingleNode() :
         vpindex(SIZE_MAX), left(SIZE_MAX), right(SIZE_MAX), radius(-INFINITY),
         sameCluster(false), maxindex(SIZE_MAX), childL(NULL), childR(NULL)  { }

   HClustVpTreeSingleNode(size_t left, size_t right) :
         vpindex(SIZE_MAX), left(left), right(right), radius(-INFINITY),
         sameCluster(false), maxindex(SIZE_MAX), childL(NULL), childR(NULL)  { }

   HClustVpTreeSingleNode(size_t vpindex, size_t left, size_t right, double radius) :
         vpindex(vpindex), left(left), right(right), radius(radius),
         sameCluster(false), maxindex(SIZE_MAX), childL(NULL), childR(NULL)  { }

   ~HClustVpTreeSingleNode() {
      if (childL) delete childL;
      if (childR) delete childR;
   }
};


class HClustVpTreeSingle : public HClustNNbasedSingle
{
protected:
   HClustVpTreeSingleNode* root;
   // bool visitAll; // for testing only

   size_t chooseNewVantagePoint(size_t left, size_t right);
   HClustVpTreeSingleNode* buildFromPoints(size_t left, size_t right, std::vector<double>& distances);

   inline void getNearestNeighborsFromMinRadiusRecursive(HClustVpTreeSingleNode* node,
      size_t index, size_t clusterIndex, double minR, std::priority_queue<double>& bestR, double& maxR, NNHeap& nnheap)
   {
      // search within (minR, maxR]
      STOPIFNOT(node != NULL);
      #ifdef GENERATE_STATS
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
         ++stats.nodeVisit;
      #endif

      if (!prefetch && node->sameCluster && clusterIndex == ds.find_set(node->left))
         return;

      if (node->vpindex == SIZE_MAX) { // leaf
         getNearestNeighborsFromMinRadiusRecursiveLeaf(node, index, clusterIndex,
            minR, bestR, maxR, nnheap);
      }
      else {
         getNearestNeighborsFromMinRadiusRecursiveNonLeaf(node, index, clusterIndex,
            minR, bestR, maxR, nnheap);
      }
   }

   void getNearestNeighborsFromMinRadiusRecursiveLeaf(HClustVpTreeSingleNode* node,
      size_t index, size_t clusterIndex, double minR, std::priority_queue<double>& bestR, double& maxR, NNHeap& nnheap);
   void getNearestNeighborsFromMinRadiusRecursiveNonLeaf(HClustVpTreeSingleNode* node,
      size_t index, size_t clusterIndex, double minR, std::priority_queue<double>& bestR, double& maxR, NNHeap& nnheap);

   virtual void getNearestNeighborsFromMinRadius(size_t index, size_t clusterIndex, double minR, NNHeap& nnheap) {
      std::priority_queue<double> bestR;
      size_t minNN = (prefetch)?opts->minNNPrefetch:opts->minNNMerge;
      for (size_t i=0; i<minNN; ++i) bestR.push(INFINITY);

      double maxR = INFINITY;

      getNearestNeighborsFromMinRadiusRecursive(root, index, clusterIndex, minR, bestR, maxR, nnheap);
   }

   void updateSameClusterFlag(HClustVpTreeSingleNode* node);

   void print(HClustVpTreeSingleNode* node);

public:

   HClustVpTreeSingle(Distance* dist, HClustOptions* opts);
   ~HClustVpTreeSingle();

   virtual void print();

}; // class

} // namespace grup


#endif
