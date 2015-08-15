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

#ifndef __HCLUST2_COMMON_H
#define __HCLUST2_COMMON_H

// ---------------------------------------------------------------------------

#define DEFAULT_MAX_LEAVES_ELEMS 4
#define DEFAULT_MAX_NN_PREFETCH 2
#define DEFAULT_VP_SELECT_SCHEME 3
#define DEFAULT_VP_SELECT_CAND 5
#define DEFAULT_VP_SELECT_TEST 12
#define DEFAULT_GNAT_DEGREE 50
#define DEFAULT_GNAT_CANDIDATES_TIMES 3
#define DEFAULT_GNAT_MIN_DEGREE 2
#define DEFAULT_GNAT_MAX_DEGREE 200
#define DEFAULT_GNAT_MAX_TIMES_DEGREE 5
#define DEFAULT_EXEMPLAR_UPDATE_METHOD 2

// ---------------------------------------------------------------------------

#include "hclust2_distance.h"
#include <queue>
#include <deque>
#include <vector>


namespace DataStructures
{
template<typename T> struct Matrix
{
   size_t nrow;
   size_t ncol;
   T* data;

   Matrix() : nrow(0), ncol(0), data(NULL) {}

   Matrix(size_t nrow, size_t ncol) :
      nrow(nrow), ncol(ncol) {
      data = new T[nrow*ncol];
   }

   ~Matrix() {
      if (data) delete [] data;
   }

   Matrix(const Matrix& m) {
      nrow = m.nrow;
      ncol = m.ncol;
      data = new T[nrow*ncol];
      for (size_t i=0; i<nrow*ncol; ++i)
         data[i] = m.data[i];
   }

   Matrix& operator=(const Matrix& m) {
      if (data) delete [] data;
      nrow = m.nrow;
      ncol = m.ncol;
      data = new T[nrow*ncol];
      for (size_t i=0; i<nrow*ncol; ++i)
         data[i] = m.data[i];
      return *this;
   }

   inline T& operator()(size_t row, size_t col) { return data[col*nrow+row]; }
};


struct HeapNeighborItem
{
   size_t index;
   double dist;

   HeapNeighborItem(size_t index, double dist) :
      index(index), dist(dist) {}

   HeapNeighborItem() :
      index(SIZE_MAX), dist(-INFINITY) {}

   inline bool operator<( const HeapNeighborItem& o ) const {
      return dist < o.dist;
   }
};


struct HeapHierarchicalItem
{
   size_t index1;
   size_t index2;
   double dist;

   HeapHierarchicalItem(size_t index1, size_t index2, double dist) :
      index1(index1), index2(index2), dist(dist) {}

   inline bool operator<( const HeapHierarchicalItem& o ) const {
      return dist >= o.dist;
   }
};


struct NNHeap {
   std::priority_queue<HeapNeighborItem> heap;
   size_t maxNNPrefetch;
// #ifdef _OPENMP
//    omp_lock_t lock;
// #endif

   NNHeap(size_t maxNNPrefetch) :
         heap(std::priority_queue<HeapNeighborItem>()),
         maxNNPrefetch(maxNNPrefetch) {
// #ifdef _OPENMP
//      omp_init_lock(&lock);
// #endif
   }

   ~NNHeap() {
// #ifdef _OPENMP
//      omp_destroy_lock(&lock);
// #endif
   }

   inline void insert(double index, double dist, double& maxR) {
// #ifdef _OPENMP
//       omp_set_lock(&lock);
// #endif
      if (heap.size() >= maxNNPrefetch && dist < maxR) {
         while (!heap.empty() && heap.top().dist == maxR) {
            heap.pop();
         }
      }
      heap.push( HeapNeighborItem(index, dist) );
      if (heap.size() >= maxNNPrefetch) maxR = heap.top().dist;
// #ifdef _OPENMP
//       omp_unset_lock(&lock);
// #endif
   }

   inline void fill(std::deque<HeapNeighborItem>& nearestNeighbors) {
      while (!heap.empty()) {
         nearestNeighbors.push_front(heap.top());
         heap.pop();
      }
   }
};


struct DistanceComparator
{
   size_t index;
   Distance* distance;

   DistanceComparator(size_t index, Distance* distance)
      : index(index), distance(distance) {}

   inline bool operator()(size_t a, size_t b) {
      return (*distance)( index, a ) < (*distance)( index, b );
   }
};


struct DistanceComparatorCached
{
   std::vector<double>* distances;

   DistanceComparatorCached(std::vector<double>* distances)
      : distances(distances) {}

   inline bool operator()(size_t a, size_t b) {
      return (*distances)[a] < (*distances)[b];
   }
};


struct IndexComparator
{
   size_t index;

   IndexComparator(size_t index)
      : index(index) {}

   inline bool operator()(size_t a) {
      return (a <= index);
   }
};

inline bool comparer_gt(size_t i, size_t j) { return (i>j); }

struct HClustTreeOptions
{
   size_t degree;           // for GNAT
   size_t candidatesTimes;  // for GNAT
   size_t minDegree;        // for GNAT
   size_t maxDegree;        // for GNAT
   size_t maxTimesDegree;   // for GNAT
   size_t maxLeavesElems;
   size_t maxNNPrefetch;    // general
   size_t vpSelectScheme;   // vp-tree and GNAT
   size_t vpSelectCand;     // for vpSelectScheme == 1
   size_t vpSelectTest;     // for vpSelectScheme == 1
   size_t exemplarUpdateMethod; // exemplar - naive(0) or not naive(1)?

   HClustTreeOptions(Rcpp::RObject control);
   Rcpp::NumericVector toR() const;
};


struct HClustTreeStats
{
   size_t nodeCount; // how many nodes are there in the tree
   size_t leafCount; // how many leaves
   size_t nodeVisit; // how many nodes were visited during NN search
   size_t nnCals;    // how many times NN search job was launched
   size_t nnCount;   // how many NNs were obtained in overall

   HClustTreeStats();
   ~HClustTreeStats();
   Rcpp::NumericVector toR() const;
};


struct SortedPoint
{
   size_t i;
   size_t j;

   SortedPoint()
      :i(0),j(0) {}

   SortedPoint(size_t _i, size_t _j) {
      if(_j < _i) {
         i = _j;
         j = _i;
      }
      else {
         i = _i;
         j = _j;
      }
   }

   inline bool operator==(const SortedPoint &other) const {
      return (i == other.i && j == other.j);
   }
};


struct Point
{
   size_t i;
   size_t j;

   Point()
      : i(0),j(0) {}

   Point(size_t _i, size_t _j) {
         i = _i;
         j = _j;
   }

   inline bool operator==(const Point &other) const {
      return (i == other.i && j == other.j);
   }
};

} // namespace DataStructures


namespace std
{

template <> struct hash<DataStructures::Point>
{
   std::size_t operator()(const DataStructures::Point& k) const {
     std::size_t seed = 0;
     boost::hash_combine(seed, k.i);
     boost::hash_combine(seed, k.j);
     return seed;
   }
};

template <> struct hash<DataStructures::SortedPoint>
{
   std::size_t operator()(const DataStructures::SortedPoint& k) const {
     std::size_t seed = 0;
     boost::hash_combine(seed, k.i);
     boost::hash_combine(seed, k.j);
     return seed;
   }
};

} // namespace std



#endif
