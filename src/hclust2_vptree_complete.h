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


#ifndef __HCLUST2_VPTREE_COMPLETE_H
#define __HCLUST2_VPTREE_COMPLETE_H
// ************************************************************************

//#define MB_IMPROVEMENT
//#define USE_BOOST_DISJOINT_SETS


// ************************************************************************
#define VANTAGE_POINT_SELECT_SCHEME 3


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

#ifdef MB_IMPROVEMENT
#include <unordered_set>
#endif // MB_IMPROVEMENT
namespace DataStructures
{

struct HClustBiVpTreeCompleteNode
{
   size_t vpindex;
   size_t left;
   size_t right;
   double radius;
   bool sameCluster;
   HClustBiVpTreeCompleteNode *ll, *lr, *rl, *rr;

   HClustBiVpTreeCompleteNode() :
      vpindex(SIZE_MAX), left(SIZE_MAX), right(SIZE_MAX), radius(-INFINITY),
      sameCluster(false), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

   HClustBiVpTreeCompleteNode(size_t left, size_t right) :
      vpindex(SIZE_MAX), left(left), right(right), radius(-INFINITY),
      sameCluster(false), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

   HClustBiVpTreeCompleteNode(size_t vpindex, double radius) :
      vpindex(vpindex), left(SIZE_MAX), right(SIZE_MAX), radius(radius),
      sameCluster(false), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

   ~HClustBiVpTreeCompleteNode() {
      if(ll) delete ll;
      if(lr) delete lr;
      if(rl) delete rl;
      if(rr) delete rr;
   }
};

class HClustBiVpTreeComplete
{
protected:

   struct HeapHierarchicalItemMax {
         size_t index1;
         size_t index2;
         double dist;
         size_t iter;

         HeapHierarchicalItemMax(size_t index1, size_t index2, double dist, size_t iter) :
            index1(index1), index2(index2), dist(dist), iter(iter) {}

         bool operator<( const HeapHierarchicalItemMax& o ) const {
            return dist >= o.dist;
         }
      };

   struct KKItem {
            double dist;
            size_t iter;

            KKItem(double dist, size_t iter) :
               dist(dist), iter(iter) {}
            KKItem() :
               dist(INFINITY), iter(0) {}
         };

   size_t maxNumberOfElementsInLeaves; // set in the constructor
   const size_t maxNearestNeighborPrefetch = 1;

   HClustBiVpTreeCompleteNode* _root;
   size_t _n;
   Distance* _distance;
   std::vector<size_t> _indices;

   std::vector<size_t> neighborsCount;
   std::vector<double> minRadiuses;
   std::vector<double> maxRadiuses;
   std::vector<bool> shouldFind;
   std::vector< deque<HeapNeighborItem> > nearestNeighbors;

   unordered_map<SortedPoint, KKItem> KK;
   unordered_map<size_t, size_t> timestamp;

   std::map<size_t,size_t> rank;
   std::map<size_t,size_t> parent;

#ifdef USE_BOOST_DISJOINT_SETS
   boost::disjoint_sets<
     associative_property_map< std::map<size_t,size_t> >,
     associative_property_map< std::map<size_t,size_t> > > ds;
#else
   PhatDisjointSets ds;
#endif
#ifdef MB_IMPROVEMENT
   unordered_map<SortedPoint, double> distClust;
   unordered_set<size_t> clusters;
   bool mbimprovement = false;
#endif


   int chooseNewVantagePoint(size_t left, size_t right);


   HClustBiVpTreeCompleteNode* buildFromPoints(size_t left, size_t right);

   void getNearestNeighborsFromMinRadiusRecursive( HClustBiVpTreeCompleteNode* node, size_t index,
      size_t clusterIndex, double minR, double& maxR,
      std::priority_queue<HeapNeighborItem>& heap );

   void print(HClustBiVpTreeCompleteNode* n);

public:
   HeapNeighborItem getNearestNeighbor(size_t index);

#ifdef DEBUG
   void printCounters();
#endif

public:

   // constructor (OK, we all know what this is, but I label it for faster in-code search)
   HClustBiVpTreeComplete(Distance* dist, size_t maxNumberOfElementsInLeaves);
   virtual ~HClustBiVpTreeComplete();
   void print();
   HeapHierarchicalItemMax calculateCluster2ClusterMaxDistance(size_t item1, size_t item2, size_t iter);
   NumericMatrix compute();


}; // class
}//namespace
#endif
