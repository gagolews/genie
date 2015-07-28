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

// #define USE_BOOST_DISJOINT_SETS
#define USE_ONEWAY_VPTREE

// ************************************************************************

#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <queue>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
// #include <fstream>
// #include <deque>
// #include <exception>
// #include <string>
// #include <boost/property_map/property_map.hpp>
// #include <boost/tuple/tuple_comparison.hpp>
// #include <algorithm>

#include "hclust2_common.h"
#include "hclust2_merge.h"
#ifdef USE_BOOST_DISJOINT_SETS
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#else
#include "disjoint_sets.h"
#endif


namespace DataStructures
{

#ifdef USE_ONEWAY_VPTREE
   #define CHILD_L  0
   #define CHILD_R  1
   #define CHILD_NUM 2
#else
   #define CHILD_LL 0
   #define CHILD_LR 1
   #define CHILD_RL 2
   #define CHILD_RR 3
   #define CHILD_NUM 4
#endif

struct HClustBiVpTreeSingleNode
{
   size_t vpindex;
   size_t left;
   size_t right;
   double radius;
   bool sameCluster;
#ifndef USE_ONEWAY_VPTREE
   HClustBiVpTreeSingleNode* children[4];
#else
   HClustBiVpTreeSingleNode* children[2];
   size_t maxindex;
#endif

   HClustBiVpTreeSingleNode() :
         vpindex(SIZE_MAX), left(SIZE_MAX), right(SIZE_MAX), radius(-INFINITY),
         sameCluster(false)  {
#ifdef USE_ONEWAY_VPTREE
      maxindex = SIZE_MAX;
      children[CHILD_L] = NULL;
      children[CHILD_R] = NULL;
#else
      children[CHILD_LL] = NULL;
      children[CHILD_LR] = NULL;
      children[CHILD_RL] = NULL;
      children[CHILD_RR] = NULL;
#endif
      }

   HClustBiVpTreeSingleNode(size_t left, size_t right) :
         vpindex(SIZE_MAX), left(left), right(right), radius(-INFINITY),
         sameCluster(false)   {
#ifdef USE_ONEWAY_VPTREE
      maxindex = SIZE_MAX;
      children[CHILD_L] = NULL;
      children[CHILD_R] = NULL;
#else
      children[CHILD_LL] = NULL;
      children[CHILD_LR] = NULL;
      children[CHILD_RL] = NULL;
      children[CHILD_RR] = NULL;
#endif
      }

   HClustBiVpTreeSingleNode(size_t vpindex, double radius) :
         vpindex(vpindex), left(SIZE_MAX), right(SIZE_MAX), radius(radius),
         sameCluster(false)   {
#ifdef USE_ONEWAY_VPTREE
      maxindex = SIZE_MAX;
      children[CHILD_L] = NULL;
      children[CHILD_R] = NULL;
#else
      children[CHILD_LL] = NULL;
      children[CHILD_LR] = NULL;
      children[CHILD_RL] = NULL;
      children[CHILD_RR] = NULL;
#endif
      }

   ~HClustBiVpTreeSingleNode() {
#ifdef USE_ONEWAY_VPTREE
      if (children[CHILD_L])  delete children[CHILD_L];
      if (children[CHILD_R])  delete children[CHILD_R];
#else
      if (children[CHILD_LL]) delete children[CHILD_LL];
      if (children[CHILD_LR]) delete children[CHILD_LR];
      if (children[CHILD_RL]) delete children[CHILD_RL];
      if (children[CHILD_RR]) delete children[CHILD_RR];
#endif
   }
};


class HClustBiVpTreeSingle
{
protected:

   HClustBiVpTreeOptions opts;

   HClustBiVpTreeSingleNode* _root;
   size_t _n;
   Distance* _distance;
   std::vector<size_t> _indices;

   std::vector<size_t> neighborsCount;
   std::vector<double> minRadiuses;
   // std::vector<double> maxRadiuses;
   std::vector<bool> shouldFind;
   std::vector< deque<HeapNeighborItem> > nearestNeighbors;

   std::map<size_t,size_t> rank;
   std::map<size_t,size_t> parent;

   HClustBiVpTreeStats stats;

#ifdef USE_BOOST_DISJOINT_SETS
   boost::disjoint_sets<
     associative_property_map< std::map<size_t,size_t> >,
     associative_property_map< std::map<size_t,size_t> > > ds;
#else
   PhatDisjointSets ds;
#endif

   int chooseNewVantagePoint(size_t left, size_t right);
   HClustBiVpTreeSingleNode* buildFromPoints(size_t left, size_t right);

   void getNearestNeighborsFromMinRadiusRecursive(HClustBiVpTreeSingleNode* node,
      size_t index, size_t clusterIndex, double minR, double& maxR,
      std::priority_queue<HeapNeighborItem>& heap);

   void print(HClustBiVpTreeSingleNode* n);


public:

   HClustBiVpTreeSingle(Distance* dist, RObject control);
   ~HClustBiVpTreeSingle();

   void print();
   NumericMatrix compute();

   HeapNeighborItem getNearestNeighbor(size_t index);

   inline const HClustBiVpTreeStats& getStats() { return stats; }
   inline const HClustBiVpTreeOptions& getOptions() { return opts; }

}; // class

} // namespace DataStructures


#endif
