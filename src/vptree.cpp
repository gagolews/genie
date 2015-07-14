#ifndef VPTREE_H_
#define VPTREE_H_
#include <iostream>
#include <Rcpp.h>
#define USE_RINTERNALS
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <queue>
#include <limits>
#include <numeric>
#include <fstream>
#include <exception>
#include <unordered_map>
#include <string>

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "boost/serialization/unordered_map.hpp"
#include "metric_trees_helpers.h"

using namespace Rcpp;
using namespace std;

template<typename T>
   class VpTree
{
   public:
   static const char* ClassName;

   VpTree(int m=16, int minM=8, size_t max_leaf_size=10, size_t vantage_point_candidates=3, size_t test_point_count=5)
      : _m(m), _minM(minM), _max_leaf_size(max_leaf_size), _vantage_point_candidates(vantage_point_candidates),
   _test_point_count(test_point_count), _root(NULL)
   {
      if(max_leaf_size < vantage_point_candidates + test_point_count)
      {
         stop("Error: max leaf size is too small");
      }
   }

   virtual ~VpTree() {
      for(size_t i=0;i<_items.size();i++)
         R_ReleaseObject(_items[i]);
      delete _root;
   }

   vector<T> getCopyItems()
   {
      return vector<T>(_items.begin(), _items.end());
   }

   Function getMetricFunction()
   {
      return _distance.distance->f;
   }

   size_t treeSize()
   {
      if(_root==NULL) return sizeof(VpTree);
      return sizeof(VpTree) + treeSize_rec(_root);
   }

   int treeHeight()
   {
      if(_root==NULL) return 0;
      return treeHeight_rec(_root);
   }

   void create( const std::vector<T>& items) {
      //stop("tworze drzewo");
      //Rcout << "tworze drzewo z " << items.size() << "elementow" << endl;
      delete _root;
      //_distance = distance;
      _items = items;
      _distance.items = &_items;
      #ifdef DEBUG
      _distance.metricCalculated = 0;
      _distance.hashmapHit = 0;
      #endif
      vector<int> indices(items.size());
      for(size_t i=0;i<indices.size();i++)
         indices[i] = i;
      //_root = buildFromPoints(indices, _minM, _minM, NULL, NULL);
      _root = buildFromPoints(indices, _m, _m, NULL, NULL);
      //Rcout << "items size = " << _items.size() << endl;
      #ifdef DEBUG
      Rcout << "metric calculated = " << _distance.metricCalculated << endl;
      Rcout << "hashmapHit = " << _distance.hashmapHit << endl;
      #endif
   }

   void setItems(const std::vector<T>& items)
   {
      _items = items;
      _distance.items = &_items;
   }

   void setDistanceFunction(const distClass& distance_arg)
   {
      _distance = distance_arg;
      _distance.items = &_items;
      //isSimilarity = _distance.isSimilarity;
   }

   void setMetricFunction(Function f)
   {
      //distClass distci;
      //distci.isSimilarity = isSimilarity;
      _distance.distance = new RFunction(f);
      //_distance = distci;
      //_distance.items = &_items;
   }

   void searchKNN( const T& target, int k, std::vector<T>* results,
                  std::vector<double>* distances, bool findItself = true)
   {
      std::priority_queue<HeapItem> heap;

      _tau = std::numeric_limits<double>::max();
      search( _root, target, true, k, heap, findItself );

      results->clear(); distances->clear();

      while( !heap.empty() ) {
         results->push_back( _items[heap.top().index] );
         distances->push_back( heap.top().dist );
         heap.pop();
      }

      std::reverse( results->begin(), results->end() );
      std::reverse( distances->begin(), distances->end() );
   }

   void searchKNNKnown( const T& target, int k, std::vector<T>* results,
                       std::vector<double>* distances, bool findItself = true)
   {
      int index = findIndex(target);
      searchKNNKnownIndex(index, k, results, distances, findItself);
   }

   virtual void searchKNNKnownIndex(int index, int k, std::vector<T>* results,
                            std::vector<double>* distances, bool findItself = true)
   {
      if(index < 0 || index >= _items.size()) stop("Index out of bounds.");
      #ifdef DEBUG
      _distance.metricCalculated = 0;
      _distance.hashmapHit = 0;
      #endif
      std::priority_queue<HeapItem> heap;

      _tau = std::numeric_limits<double>::max();
      search( _root, index, true, k, heap, findItself );

      results->clear(); distances->clear();

      while( !heap.empty() ) {
         results->push_back( _items[heap.top().index] );
         distances->push_back( heap.top().dist );
         heap.pop();
      }

      std::reverse( results->begin(), results->end() );
      std::reverse( distances->begin(), distances->end() );
      #ifdef DEBUG
      Rcout << "metric calculated = " << _distance.metricCalculated << endl;
      Rcout << "hashmapHit = " << _distance.hashmapHit << endl;
      #endif
   }

   void searchRadiusKnownIndex(int index, double tau, std::vector<T>* results,
                               std::vector<double>* distances, bool findItself = true)
   {
      if(index < 0 || index >= _items.size()) stop("Index out of bounds.");
      std::priority_queue<HeapItem> heap;

      _tau = tau;
      search( _root, index, false, -1, heap, findItself );

      results->clear(); distances->clear();

      while( !heap.empty() ) {
         results->push_back( _items[heap.top().index] );
         distances->push_back( heap.top().dist );
         heap.pop();
      }

      std::reverse( results->begin(), results->end() );
      std::reverse( distances->begin(), distances->end() );
   }

   void searchRadiusKnown( const T& target, double tau, std::vector<T>* results,
                          std::vector<double>* distances, bool findItself = true)
   {
      int index = findIndex(target);
      searchRadiusKnownIndex(index, tau, results, distances, findItself);
   }

   void searchRadius( const T& target, double tau, std::vector<T>* results,
                     std::vector<double>* distances, bool findItself = true)
   {
      std::priority_queue<HeapItem> heap;

      _tau = tau;
      search( _root, target, false, -1, heap, findItself);

      results->clear(); distances->clear();

      while( !heap.empty() ) {
         results->push_back( _items[heap.top().index] );
         distances->push_back( heap.top().dist );
         heap.pop();
      }

      std::reverse( results->begin(), results->end() );
      std::reverse( distances->begin(), distances->end() );
   }

   void insert(const T& target)
   {
      _items.push_back(target);
      insert_rec(_root, target);
   }

   #ifdef DEBUG
   void printCounters()
   {
      _distance.printCounters();
   }
   #endif

   protected:
   distClass _distance;
   std::vector<T> _items;
   double _tau;

   // The following condition must held:
   // MAX_LEAF_SIZE >= VANTAGE_POINT_CANDIDATES + TEST_POINT_COUNT
   int _m; //max number of children
   int _minM; //min number of children
   size_t _max_leaf_size;
   size_t _vantage_point_candidates;
   size_t _test_point_count;

   friend class boost::serialization::access;
   template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
   {
      //ar & _items;
      //ar & _tau;
      ar & _m;
      ar & _minM;
      ar & _max_leaf_size;
      ar & _vantage_point_candidates;
      ar & _test_point_count;
      ar & _distance;
      ar & _root;
   }

   struct Node
   {
      int vpindex;
      //Node* left;
      //Node* right;
      bool isLeaf;
      vector<double> radiuses;
      vector<int> points;
      vector<Node*> children;
      int childCount;
      Node* parent;

      friend class boost::serialization::access;
      template<class Archive>
         void serialize(Archive & ar, const unsigned int version)
      {
         ar & vpindex;
         //for(int i=0; i<children.size(); i++)
         //	ar & children[i];
         ar & children;
         ar & parent;
         ar & points;
         ar & radiuses;
         ar & isLeaf;
         ar & childCount;
      }

      Node() :
      vpindex(0), isLeaf(false), childCount(0) {}

      ~Node() {
         for(int i=0;i<childCount;i++)
         {
            if(children[i] != NULL)
               delete children[i];
         }
      }
   }* _root;

   struct HeapItem {
      HeapItem( int index, double dist) :
      index(index), dist(dist) {}
      int index;
      double dist;
      bool operator<( const HeapItem& o ) const {
         return dist < o.dist;
      }
   };

   struct DistanceComparator
   {
      vector<T>* items;
      int index;
      distClass& distance;

      //_items, indices[0], _distance
      DistanceComparator(vector<T>* items, int index, distClass& distance ) : items(items), index(index), distance(distance) {}
      bool operator()(int a, int b) {
         return distance( index, a ) < distance( index, b );
      }
   };

   int findIndex(const T& target)
   {
      //      int a = 43;
      //NumericVector nvTarget(target);
      for(size_t i = 0; i<_items.size(); i++)
      {
         //if(Rcpp::all(_items[i] == target))
         //NumericVector nvTarget(target);
         if(R_compute_identical(_items[i], target, 16))
            return i;
      }
      stop("There is no such element in the tree.");
      return -1; // fake retval, avoid warning
   }

   void splitLeaf(Node* node, const T& target)
   {
      Node* parent = node->parent;
      if(parent == NULL) //strange case, when we have only one leaf in whole tree, very small number of items
      {
         //Rcout << "splitLeaf - tree grows in height, items size = " << _items.size() << endl;
         delete _root;
         vector<int> indices(_items.size());
         for(size_t i=0;i<indices.size();i++)
            indices[i] = i;
         _root = buildFromPoints(indices, _minM, _minM, NULL);
         return;
      }
      //Rcout << "splitLeaf - przebudowuje" << endl;
      vector<int> allitems;
      allitems.push_back(_items.size() - 1);

      for(int i=0;i<parent->childCount;i++)
      {
         vector<int> items = parent->children[i]->points;
         allitems.insert(allitems.end(), items.begin(), items.end());
      }

      buildFromPoints(allitems, parent->childCount+1, _minM, parent->parent, parent);
   }

   void splitNonLeaf(Node* node, const T& target)
   {
      Node* parent = node->parent;
      while(parent != NULL && parent->childCount >= _m)
         parent = parent->parent;

      if(parent == NULL) // If root is full, then tree grows in height
      {
         //Rcout << "splitNonLeaf - tree grows in height, items size = " << _items.size() << endl;
         delete _root;
         vector<int> indices(_items.size());
         for(size_t i=0;i<indices.size();i++)
            indices[i] = i;
         _root = buildFromPoints(indices, _minM, _minM, NULL, NULL);
         return;
      }
      //Rcout << "splitNonLeaf - przebudowuje" << endl;
      vector<int> allitems = getAllItemsFromNode(parent);
      allitems.push_back(_items.size() - 1);
      buildFromPoints(allitems, parent->childCount+1, _minM, parent->parent, parent);
   }

   vector<int> getAllItemsFromNode(Node* node)
   {
      if(node->isLeaf)
         return node->points;

      vector<int> ret;
      for(int i=0;i<node->childCount;i++)
      {
         vector<int> items = getAllItemsFromNode(node->children[i]);
         ret.insert(ret.end(), items.begin(), items.end());
      }
      return ret;
   }

   void insert_rec(Node* node, const T& target)
   {
      if ( node == NULL ) return;
      if(node->isLeaf)
      {
         if(node->points.size() < _max_leaf_size)
         {
            //Rcout << "dodaje normalnie, itemow w drzewie: " << _items.size() << endl;
            node->points.push_back(_items.size() - 1);
            //Rcout << _items[_items.size() - 1][0] << ", " << _items[_items.size() - 1][1] << endl;
            /*for(int i=0;i<node->points.size(); i++)
    				Rcout << node->points[i] << " ";
    			Rcout << endl;*/
         }
         else
         {
            if(node->parent->childCount < _m)
            {
               //Rcout << "splitLeaf" << endl;
               splitLeaf(node, target);
            }
            else
            {
               //Rcout << "splitNonLeaf" << endl;
               splitNonLeaf(node, target);
            }
         }
      }
      else
      {
         double dist = _distance(_items[node->vpindex], target);
         for(int i=0;i<node->childCount-1;i++)
         {
            if ( dist < node->radiuses[i] ) {
               insert_rec(node->children[i], target);
               return;
            }
         }
         insert_rec( node->children[node->childCount-1], target);
      }
   }

   Node* buildFromPoints(vector<int> indices, int childCountThis, int childCountChildren, Node* parent, Node* node=NULL )
   {
      if(node == NULL)
         node = new Node();
      else
      {
         for(int i=0;i<node->childCount;i++)
            delete node->children[i];
      }
      node->parent = parent;

      if ( indices.size() < _max_leaf_size  ) {
         node->isLeaf=true;
         node->points = indices;
         return node;
      }

      node->isLeaf=false;
      node->childCount = childCountThis;

      // choose an arbitrary point and move it to the start
      int vpi = chooseNewVantagePoint(indices); //(int)((double)rand() / RAND_MAX * (upper - lower - 1) ) + lower;
      std::swap( indices[0], indices[vpi] );
      std::sort(indices.begin() + 1, indices.end(), DistanceComparator( &_items, indices[0], _distance ));

      node->radiuses = vector<double>(childCountThis-1);
      node->children = vector<Node*>(childCountThis);
      for(int i=0;i<childCountThis;i++)
      {
         int itemsPerNode = indices.size() / childCountThis;
         if(i < childCountThis - 1)
         {
            node->children[i] = buildFromPoints(
               vector<int>(indices.begin() + i*itemsPerNode, indices.begin() + (i+1)*itemsPerNode),
               childCountChildren, childCountChildren,
               node);
            node->radiuses[i] = (_distance( _items[indices[0]], _items[indices[0 + (i+1)*itemsPerNode - 1] ])
                                 + _distance( _items[indices[0]], _items[indices[0+ (i+1)*itemsPerNode] ])) / 2.0;
         }
         else
         {
            node->children[i] = buildFromPoints(
               vector<int>(indices.begin() + (childCountThis - 1)*itemsPerNode, indices.end()),
               childCountChildren, childCountChildren,
               node);
         }
      }
      node->vpindex = indices[0];
      return node;
   }

   int chooseNewVantagePoint(vector<int> indices)
   {
      vector<int> candidates(0);
      vector<int> testPoints(0);

      for (size_t i = 0; i < 0 + _vantage_point_candidates; ++i) {
         int basePointIndex = i + (int) (rand() % (indices.size() - i));
         candidates.push_back(indices[basePointIndex]);
         std::swap(indices[i], indices[basePointIndex] );
      }
      //Rcout << "wybralem kandydatow" << endl;
      for (size_t i = 0 + _vantage_point_candidates; i < 0 + _vantage_point_candidates + _test_point_count; ++i) {
         int testPointIndex = i + (int) (rand() % (indices.size() - i));
         testPoints.push_back(indices[testPointIndex]);
         std::swap(indices[i], indices[testPointIndex] );
      }
      //Rcout << "wybralem testerow" << endl;
      /*if(candidates.size() != VANTAGE_POINT_CANDIDATES)
            Rcout << "candidates.size() " << candidates.size() << endl;
        if(testPoints.size() != TEST_POINT_COUNT)
            Rcout << "testPoints.size() " << testPoints.size() << endl;*/
      double bestBasePointSigma = 0;
      int bestIndex = 0;
      for(size_t i=0;i<candidates.size();i++) {
         vector<double> distances(_test_point_count);
         for (size_t j = 0; j < _test_point_count; ++j) {
            distances[j] = _distance( candidates[i], testPoints[j] );
         }
         //Rcout << "przed liczeniem sigmy" << endl;
         double sigma = sigmaSquare(distances);
         //Rcout << "policzylem sigme" << endl;
         if (sigma > bestBasePointSigma) {
            bestBasePointSigma = sigma;
            bestIndex = i;
         }
      }

      return bestIndex;
   }

   static double sigmaSquare(const vector<double>& values) {
      double sum = std::accumulate(values.begin(),values.end(),0);

      double avg = sum / values.size();
      double sigmaSq = 0;

      for (size_t i=0;i<values.size();i++) {
         double value = values[i];
         double dev = value - avg;
         sigmaSq += dev * dev;
      }

      return sigmaSq;
   }

   size_t calculateNodeSize(Node* node)
   {
      return node->radiuses.size()*sizeof(double)
         + node->points.size()*sizeof(int)
         + node->children.size()*sizeof(Node*)
         + sizeof(Node);
   }

   size_t treeSize_rec(Node* node)
   {
      size_t size = calculateNodeSize(node);
      for(int i=0;i<node->childCount;i++)
      {
         size += treeSize_rec(node->children[i]);
      }
      return size;
   }

   int treeHeight_rec(Node* node)
   {
      int maxH = 0;
      for(int i=0;i<node->childCount;i++)
      {
         maxH = max(treeHeight_rec(node->children[i]), maxH);
      }
      return maxH+1;
   }

   void search( Node* node, const T& target, bool isKNN, int k,
               std::priority_queue<HeapItem>& heap, bool findItself)
   {
      if ( node == NULL ) return;


      //printf("dist=%g tau=%gn", dist, _tau );
      if(node->isLeaf)
      {
         for(size_t i=0;i<node->points.size();i++)
         {
            double dist2 = _distance( _items[node->points[i]], target );
            if(dist2 < 1e-6)
            {
                if(!findItself && R_compute_identical(_items[node->points[i]], target, 16))
                   continue;
            }
            if ( (dist2 < _tau && isKNN) || (dist2 <= _tau && !isKNN) ) {
               if ( heap.size() == (size_t)k && isKNN) heap.pop();
               heap.push( HeapItem(node->points[i], dist2) );
               if ( heap.size() ==(size_t) k && isKNN)
               {
                  _tau = heap.top().dist;
                  //Rcout << "current tau=" << _tau << endl;
               }
            }
         }
      }
      else
      {
         double dist = _distance( _items[node->vpindex], target );
         vector<bool> visited(node->childCount, false);

         for(int i=0;i<node->childCount-1;i++)
         {
            //if ( dist < node->radiuses[i] ) {
            if ( dist - _tau <= node->radiuses[i] && !visited[i]) {
               search( node->children[i], target, isKNN, k, heap, findItself );
               visited[i]=true;
            }
            if ( dist + _tau >= node->radiuses[i] && !visited[i+1]) {
               search( node->children[i+1], target, isKNN, k, heap, findItself );
               visited[i+1]=true;
            }
            //}
         }
         /*if ( dist + _tau >= node->radiuses[node->childCount-2] ) {
				search( node->children[node->childCount-1], target, k, heap );
			}*/
      }
   }

   virtual void search( Node* node, int index, bool isKNN, int k,
               std::priority_queue<HeapItem>& heap, bool findItself )
   {
      if ( node == NULL ) return;

      //printf("dist=%g tau=%gn", dist, _tau );
      if(node->isLeaf)
      {
         for(size_t i=0;i<node->points.size();i++)
         {
            double dist2 = _distance(node->points[i], index );
            if(dist2 < 1e-6)
            {
                if(!findItself && R_compute_identical(_items[node->points[i]], _items[index], 16))
                   continue;
            }
            if ( (dist2 < _tau && isKNN) || (dist2 <= _tau && !isKNN) ) {
               if ( heap.size() ==(size_t) k && isKNN) heap.pop();
               heap.push( HeapItem(node->points[i], dist2) );
               if ( heap.size() == (size_t) k && isKNN)
               {
                  _tau = heap.top().dist;
                  //Rcout << "current tau=" << _tau << endl;
               }
            }
         }
      }
      else
      {
         double dist = _distance(node->vpindex, index);
         vector<bool> visited(node->childCount, false);

         for(int i=0;i<node->childCount-1;i++)
         {
            //if ( dist < node->radiuses[i] ) {
            if ( dist - _tau <= node->radiuses[i] && !visited[i]) {
               search( node->children[i], index, isKNN, k, heap, findItself );
               visited[i]=true;
            }
            if ( dist + _tau >= node->radiuses[i] && !visited[i+1]) {
               search( node->children[i+1], index, isKNN, k, heap, findItself );
               visited[i+1]=true;
            }
            //}
         }
         /*if ( dist + _tau >= node->radiuses[node->childCount-2] ) {
				search( node->children[node->childCount-1], target, k, heap );
			}*/
      }
   }
};
/*
template <>
   void vptree<RObject>::findIndex(const RObject& target) // specialize only one member
   {
      Rcout << "specialized" << endl;
      for(int i = 0; i<_items.size(); i++)
      {
         if(Rcpp::all(_items[i] == target))
            //if(_items[i] == target)
            return i;
      }
      stop("There is no such element in the tree.");
   }
*/
template<typename T>
   const char* VpTree<T>::ClassName = "VpTree";

template<typename T>
   void checkIsVpTreeClass(XPtr< VpTree<T> >& _tree)
{
   if(strcmp(_tree.attr("class"), VpTree<T>::ClassName))
      stop("not a VpTree object");
}


#endif /* VPTREE_H_ */
