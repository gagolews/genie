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

#ifndef __HCLUST2_VPTREE_SINGLE_H
#define __HCLUST2_VPTREE_SINGLE_H



// ************************************************************************

#include "hclust2_nnbased_single.h"


namespace DataStructures
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

   HClustVpTreeSingleNode* _root;
   std::vector<double> distances;

   size_t chooseNewVantagePoint(size_t left, size_t right);
   HClustVpTreeSingleNode* buildFromPoints(size_t left, size_t right);

   void getNearestNeighborsFromMinRadiusRecursive(HClustVpTreeSingleNode* node,
      size_t index, size_t clusterIndex, double minR, double& maxR, NNHeap& nnheap);

   virtual void getNearestNeighborsFromMinRadius(size_t index, size_t clusterIndex, double minR, NNHeap& nnheap) {
      double maxR = INFINITY;
      getNearestNeighborsFromMinRadiusRecursive(_root, index, clusterIndex, minR, maxR, nnheap);
   }

   void print(HClustVpTreeSingleNode* n);

public:

   HClustVpTreeSingle(Distance* dist, RObject control);
   ~HClustVpTreeSingle();

   virtual void print();

}; // class

} // namespace DataStructures


#endif
