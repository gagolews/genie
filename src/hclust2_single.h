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

#ifndef __HCLUST2_SINGLE_H
#define __HCLUST2_SINGLE_H



// ************************************************************************

// #define DEFAULT_MAX_NUM_ELEMS_LEAVES 2
#define DEFAULT_MAX_NN_PREFETCH 2
#define VANTAGE_POINT_SELECT_SCHEME 3
#define VANTAGE_POINT_SELECT_SCHEME_1_NUMCANDIDATES 5
#define VANTAGE_POINT_SELECT_SCHEME_1_NUMTEST 12
// #define MB_IMPROVEMENT
// #define USE_BOOST_DISJOINT_SETS
#define NN_COUNTERS


// ************************************************************************

#if VERBOSE > 7 && !defined(NN_COUNTERS)
#define NN_COUNTERS
#endif


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


#include "hclust2_distance.h"
#include "hclust2_merge.h"
#ifdef USE_BOOST_DISJOINT_SETS
#include <boost/pending/disjoint_sets.hpp>
#else
#include "disjoint_sets.h"
#endif

#ifdef MB_IMPROVEMENT
#include <unordered_set>
#endif // MB_IMPROVEMENT



namespace DataStructures{

class HClustSingleBiVpTree
{
protected:

   struct HeapNeighborItem {
      size_t index;
      double dist;

      HeapNeighborItem(size_t index, double dist) :
         index(index), dist(dist) {}

      HeapNeighborItem() :
         index(SIZE_MAX), dist(-INFINITY) {}

      bool operator<( const HeapNeighborItem& o ) const {
         return dist < o.dist;
      }
   };

   struct HeapHierarchicalItem {
      size_t index1;
      size_t index2;
      double dist;

      HeapHierarchicalItem(size_t index1, size_t index2, double dist) :
         index1(index1), index2(index2), dist(dist) {}

      bool operator<( const HeapHierarchicalItem& o ) const {
         return dist >= o.dist;
      }
   };

   struct Node
   {
      size_t vpindex;
      size_t left;
      size_t right;
      double radius;
      bool sameCluster;
      Node *ll, *lr, *rl, *rr;

      Node() :
         vpindex(SIZE_MAX), left(SIZE_MAX), right(SIZE_MAX), radius(-INFINITY),
         sameCluster(false), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

      Node(size_t left, size_t right) :
         vpindex(SIZE_MAX), left(left), right(right), radius(-INFINITY),
         sameCluster(false), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

      Node(size_t vpindex, double radius) :
         vpindex(vpindex), left(SIZE_MAX), right(SIZE_MAX), radius(radius),
         sameCluster(false), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

      ~Node() {
         if(ll) delete ll;
         if(lr) delete lr;
         if(rl) delete rl;
         if(rr) delete rr;
      }
   };

   struct DistanceComparator
   {
      size_t index;
      Distance* distance;

      DistanceComparator(size_t index, Distance* distance )
         : index(index), distance(distance) {}

      bool operator()(size_t a, size_t b) {
         return (*distance)( index, a ) < (*distance)( index, b );
      }
   };

   struct IndexComparator
   {
      size_t index;

      IndexComparator(size_t index)
         : index(index) {}

      bool operator()(size_t a) {
         return a <= index;
      }
   };

   size_t maxNumberOfElementsInLeaves; // set in the constructor
   const size_t maxNearestNeighborPrefetch = DEFAULT_MAX_NN_PREFETCH;

   Node* _root;
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

#ifdef USE_BOOST_DISJOINT_SETS
   boost::disjoint_sets<
     associative_property_map< std::map<size_t,size_t> >,
     associative_property_map< std::map<size_t,size_t> > > ds;
#else
   PhatDisjointSets ds;
#endif

#ifdef MB_IMPROVEMENT
   std::unordered_map<SortedPoint, double> distClust;
   std::unordered_set<size_t> clusters;
   bool mbimprovement = false;
#endif  // MB_IMPROVEMENT


   int chooseNewVantagePoint(size_t left, size_t right);
   Node* buildFromPoints(size_t left, size_t right);

   void getNearestNeighborsFromMinRadiusRecursive( Node* node, size_t index,
      size_t clusterIndex, double minR, double& maxR,
      std::priority_queue<HeapNeighborItem>& heap );

   void print(Node* n);


public:

   HClustSingleBiVpTree(Distance* dist, size_t maxNumberOfElementsInLeaves);
   ~HClustSingleBiVpTree();

   void print();
   NumericMatrix compute();

   HeapNeighborItem getNearestNeighbor(size_t index);

}; // class

} // namespace DataStructures


#endif
