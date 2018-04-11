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

#ifndef __HCLUST2_COMMON_H
#define __HCLUST2_COMMON_H

#include "defs.h"
#include "hclust2_distance.h"
#include "disjoint_sets.h"
#include <queue>
#include <deque>
#include <vector>
#include <queue>
#include <list>


namespace grup
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


struct HeapNeighborItemFromSmallestComparator
{
   inline bool operator()(const HeapNeighborItem& o1, const HeapNeighborItem& o2) const {
      return o1.dist >= o2.dist;
   }

};


typedef std::priority_queue<HeapNeighborItem,
   std::vector<HeapNeighborItem>, HeapNeighborItemFromSmallestComparator>
      priority_queue_HeapNeighborItem_FromSmallest;


struct HeapHierarchicalItem
{
   size_t index1;
   size_t index2;
   double dist;

   HeapHierarchicalItem() :
      index1(SIZE_MAX), index2(SIZE_MAX), dist(INFINITY) {}

   HeapHierarchicalItem(size_t index1, size_t index2, double dist) :
      index1(index1), index2(index2), dist(dist) {}

   inline bool operator<( const HeapHierarchicalItem& o ) const {
      return dist >= o.dist || (dist == o.dist && index2 > o.index2); // SIZE_MAX index2 at the end of a series
   }
};


// class HclustPriorityQueue {
// private:
//    struct BSTNode {
//       BSTNode* left;
//       BSTNode* right;
//       HeapHierarchicalItem elem;
//    };
//
//    BSTNode* root;
//
//    void rotateLeft(BSTNode** root) {
//       STOPIFNOT(*root)
//       if (!(*root)->left) return;
//       BSTNode* oldroot = *root;
//       *root = oldroot->left;
//       oldroot->left = (*root)->right;
//       (*root)->right = oldroot;
//    }
//
//    void rotateRight(BSTNode** root) {
//       STOPIFNOT(*root)
//       if (!(*root)->right) return;
//       BSTNode* oldroot = *root;
//       *root = oldroot->right;
//       oldroot->right = (*root)->left;
//       (*root)->left = oldroot;
//    }
//
//    void deleteSubTree(BSTNode** root) {
//       if (!*root) return;
//       deleteSubTree(&(*root)->left);
//       deleteSubTree(&(*root)->right);
//       delete *root;
//       *root = NULL;
//    }
//
//    void deleteLeftmost(BSTNode** root) {
//       STOPIFNOT(*root)
//       if ((*root)->left) {
//          deleteLeftmost(&(*root)->left);
//       }
//       else {
//          BSTNode* delme = *root;
//          *root = delme->right;
//          delete delme;
//       }
//    }
//
//    void insert(BSTNode** root, const HeapHierarchicalItem& data) {
//       if (*root) {
//          if (data.dist < (*root)->elem.dist)
//             insert(&(*root)->left, data);
//          else
//             insert(&(*root)->right, data);
//
//          double u = unif_rand();
//          if (u < 0.33) rotateLeft(root);
//          else if (u < 0.67) rotateRight(root);
//       }
//       else {
//          *root = new BSTNode;
//          (*root)->left = NULL;
//          (*root)->right = NULL;
//          (*root)->elem = data;
//       }
//    }
//
//
// public:
//    HclustPriorityQueue(std::size_t) { root = NULL; }
//    const bool empty() const { return root == NULL; }
//    ~HclustPriorityQueue() { deleteSubTree(&root); }
//
//    const HeapHierarchicalItem& top() {
//       STOPIFNOT(root)
//       while (root->left)
//          rotateLeft(&root);
//       return root->elem;
//    }
//
//    void pop() {
//       STOPIFNOT(root)
//       deleteLeftmost(&root);
//    }
//
//    void push(const HeapHierarchicalItem& data) {
//       insert(&root, data);
//    }
//
// };

// class HclustPriorityQueue {
// private:
//
//    std::vector<std::size_t> left;
//    std::vector<std::size_t> right;
//    std::vector<std::size_t> parent;
//    std::vector<HeapHierarchicalItem> elem;
//    std::vector<std::size_t> free;
//    std::size_t occupied;
//    std::size_t root;
//    std::size_t best;
//
//    double check_sorted;
//
//    void print(std::size_t cur, std::size_t h) {
//       if (cur == SIZE_MAX) return;
//       print(left[cur], h+1);
//       std::cerr << elem[cur].dist << "(" << h << "), ";
//       print(right[cur], h+1);
//    }
//
//    void checkSorted(std::size_t cur) {
//       if (cur == SIZE_MAX) return;
//       checkSorted(right[cur]);
//       STOPIFNOT(check_sorted >= elem[cur].dist)
//       check_sorted = elem[cur].dist;
//       checkSorted(left[cur]);
//    }
//
// public:
//
//    HclustPriorityQueue(std::size_t n) :
//         left(n), right(n), parent(n), elem(n), free(n)
//    {
//       root = SIZE_MAX;
//       occupied = 0;
//       best = SIZE_MAX;
//       for (std::size_t i=0; i<n; ++i)
//          free[i] = i;
//    }
//
//    ~HclustPriorityQueue() {  }
//
//    void print() { print(root, 1); }
//    void checkSorted() { check_sorted = INFINITY; checkSorted(root); STOPIFNOT(check_sorted == elem[best].dist) }
//
//    void pop() {
//       STOPIFNOT(best != SIZE_MAX)
//       STOPIFNOT(left[best] == SIZE_MAX) // best has no left leaves
//       free[--occupied] = best;
//       STOPIFNOT(occupied >= 0)
//       if (parent[best] == SIZE_MAX) {
//          // it's a root
//          best = root = right[best];
//          parent[root] = SIZE_MAX;
//          if (best != SIZE_MAX) {
//             while (left[best] != SIZE_MAX)
//                best = left[best];
//          }
//       }
//       else { // parent[best] != SIZE_MAX
//          STOPIFNOT(left[parent[best]] == best)
//          STOPIFNOT(elem[parent[best]].dist >= elem[best].dist)
//          STOPIFNOT(elem[root].dist >= elem[best].dist)
//          if (right[best] == SIZE_MAX) {
//             left[parent[best]] = SIZE_MAX;
//             best = parent[best];
//          }
//          else { // right[best] != SIZE_MAX
//             left[parent[best]] = right[best];
//             parent[right[best]] = parent[best];
//             best = right[best];
//             while (left[best] != SIZE_MAX)
//                best = left[best];
//          }
//       }
//    }
//
//    inline const HeapHierarchicalItem& top() {
//       STOPIFNOT(best != SIZE_MAX)
//       return elem[best];
//    }
//
//    void push(const HeapHierarchicalItem& data) {
//       STOPIFNOT(occupied+1 < free.size())
//       if (root == SIZE_MAX) {
//          STOPIFNOT(occupied == 0)
//          root = best = free[occupied++];
//          right[root] = left[root] = parent[root] = SIZE_MAX;
//          elem[root] = data;
//          return;
//       }
//
//       STOPIFNOT(best != SIZE_MAX)
//       if (data.dist < elem[best].dist) {
//          STOPIFNOT(left[best] == SIZE_MAX)
//          left[best] = free[occupied++];
//          parent[left[best]] = best;
//          best = left[best];
//          right[best] = left[best] = SIZE_MAX;
//          elem[best] = data;
//          return;
//       }
//
//
//       std::size_t start = root;
//       while (true) {
//          if (data.dist < elem[start].dist) {
//             if (left[start] == SIZE_MAX) {
//                left[start] = free[occupied++];
//                parent[left[start]] = start;
//                left[left[start]] = right[left[start]] = SIZE_MAX;
//                elem[left[start]] = data;
//                return;
//             }
//             else {
//                start = left[start];
//             }
//          }
//          else {
//             if (right[start] == SIZE_MAX) {
//                right[start] = free[occupied++];
//                parent[right[start]] = start;
//                left[right[start]] = right[right[start]] = SIZE_MAX;
//                elem[right[start]] = data;
//                return;
//             }
//             else {
//                start = right[start];
//             }
//          }
//       }
//    }
//
//    inline bool empty() const { return root == SIZE_MAX; }
// };


class HclustPriorityQueue
{
   size_t n;
   size_t ncur;
   std::vector<HeapHierarchicalItem> items;
   bool heapMade;

public:
   HclustPriorityQueue(size_t n=0) :
      n(n), ncur(0), items(n), heapMade(false) { }

   const HeapHierarchicalItem& top() {
      if (!heapMade) {
         std::make_heap(items.begin(), items.begin()+ncur);
         heapMade = true;
      }
      return items[0];
   }

   void pop() {
      if (!heapMade) {
         std::make_heap(items.begin(), items.begin()+ncur);
         heapMade = true;
      }
      std::pop_heap(items.begin(), items.begin()+ncur);
      --ncur;
      STOPIFNOT(ncur >= 0);
   }

   void push(const HeapHierarchicalItem& item) {
      STOPIFNOT(ncur < n);
      items[ncur++] = item;
      if (heapMade) {
         std::push_heap(items.begin(), items.begin()+ncur);
      }
   }

   bool empty() const {
      return (ncur == 0);
   }

   void reset() { heapMade = false; }

};



struct HClustOptions
{
//    size_t degree;           // for GNAT
//    size_t candidatesTimes;  // for GNAT
//    size_t minDegree;        // for GNAT
//    size_t maxDegree;        // for GNAT
//    size_t maxTimesDegree;   // for GNAT
   size_t maxLeavesElems;   //
   size_t maxNNPrefetch;    //
   size_t maxNNMerge;       //
   size_t minNNPrefetch;    //
   size_t minNNMerge;       //
   // std::string exemplar;      //
   bool useVpTree;
   bool useMST;
   size_t vpSelectScheme;   // vp-tree and GNAT
   size_t vpSelectCand;     // for vpSelectScheme == 1
   size_t vpSelectTest;     // for vpSelectScheme == 1
   size_t nodesVisitedLimit;// for single approx
   double thresholdGini;    // for single approx
   // size_t exemplarUpdateMethod; // exemplar - naive(0) or not naive(1)?
   // size_t maxExemplarLeavesElems; //for exemplars biggers numbers are needed I think
   // bool isCurseOfDimensionality;

   HClustOptions(Rcpp::RObject control);
   Rcpp::NumericVector toR() const;
};


struct HClustStats
{
   size_t nodeCount; // how many nodes are there in the tree
   size_t leafCount; // how many leaves
   size_t nodeVisit; // how many nodes were visited during NN search
   size_t nnCals;    // how many times NN search job was launched
   size_t nnCount;   // how many NNs were obtained in overall
   size_t medoidOldNew; //..how many times it was successful
   size_t medoidUpdateCount; // how many times we calculate d_old and d_new..

   HClustStats();
   ~HClustStats();
   Rcpp::NumericVector toR() const;
};





struct NNHeap {
   std::priority_queue< HeapNeighborItem > heap;
   static HClustOptions* opts;
   size_t exemplarsCount;
// #ifdef _OPENMP
//    omp_lock_t lock;
// #endif
   NNHeap() :
         heap(),
         exemplarsCount(0) {
// #ifdef _OPENMP
//      omp_init_lock(&lock);
// #endif
   }

   static void setOptions(HClustOptions* newopts) {
      opts = newopts;
   }

   ~NNHeap() {
// #ifdef _OPENMP
//      omp_destroy_lock(&lock);
// #endif

   }

   inline bool empty()
   {
      return heap.empty();
   }

   inline const HeapNeighborItem& top()
   {
      return heap.top();
   }

   inline const size_t size()
   {
      return heap.size();
   }

   inline void pop()
   {
      heap.pop();
   }

   inline void push(const HeapNeighborItem& elem)
   {
      heap.push(elem);
   }

   inline void insert(double index, double dist, double& maxR) {
      STOPIFNOT(NNHeap::opts != NULL)
// #ifdef _OPENMP
//       omp_set_lock(&lock);
// #endif
      if (heap.size() >= opts->maxNNPrefetch && dist < maxR) {
         while (!heap.empty() && heap.top().dist == maxR) {
            heap.pop();
         }
      }
      heap.push( HeapNeighborItem(index, dist) );
      if (heap.size() >= opts->maxNNPrefetch) maxR = heap.top().dist;
// #ifdef _OPENMP
//       omp_unset_lock(&lock);
// #endif
   }

   inline void insertExemplars(double index, double dist, double& maxR, DisjointSets& ds, bool isExemplar) {
   // #ifdef _OPENMP
   //       omp_set_lock(&lock);
   // #endif
      heap.push( HeapNeighborItem(index, dist) );
      if(isExemplar)
      {
         exemplarsCount++;
      }
      std::list<HeapNeighborItem> toRemove;
      size_t toRemoveExemplarsCount=0;

      if (heap.size() >= opts->maxNNPrefetch+1 && dist < maxR) {
         while (!heap.empty() && heap.top().dist == maxR) {
            toRemove.push_back(heap.top());
            if(heap.top().index == ds.find_set(heap.top().index))
            {
               toRemoveExemplarsCount++;
            }
            heap.pop();
         }
      }

      if(toRemoveExemplarsCount == exemplarsCount && exemplarsCount > 0)
      {
         for(auto it = toRemove.begin(); it != toRemove.end(); ++it)
            heap.push(*it);
      }
      else
      {
         exemplarsCount -= toRemoveExemplarsCount;
      }

      if (heap.size() >= opts->maxNNPrefetch && exemplarsCount > 0) maxR = heap.top().dist;
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

   inline void fill(std::list<HeapNeighborItem>& nearestNeighbors) {
      while (!heap.empty()) {
         nearestNeighbors.push_front(heap.top());
         heap.pop();
      }
   }

   inline void fill(std::priority_queue<HeapNeighborItem, std::vector<HeapNeighborItem>, HeapNeighborItemFromSmallestComparator>& nearestNeighbors)
   {
      while (!heap.empty()) {
         nearestNeighbors.push(heap.top());
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
inline bool comparer_gt(double i, double j) { return (i>j); }

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

} // namespace grup



// #include <boost/functional/hash.hpp>

namespace std
{

template <> struct hash<grup::Point>
{
   /*
    * template<typename T> void hash_combine(size_t & seed, T const& v);
    * seed ^= hash_value(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    */
   std::size_t operator()(const grup::Point& k) const {
     std::size_t seed = 0;

     seed ^= (size_t)(k.i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
     seed ^= (size_t)(k.j) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

     // boost::hash_combine(seed, k.i);
     // boost::hash_combine(seed, k.j);
     return seed;
   }
};

template <> struct hash<grup::SortedPoint>
{
   std::size_t operator()(const grup::SortedPoint& k) const {
     std::size_t seed = 0;
     seed ^= (size_t)(k.i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
     seed ^= (size_t)(k.j) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
     // boost::hash_combine(seed, k.i);
     // boost::hash_combine(seed, k.j);
     return seed;
   }
};

} // namespace std



#endif
