#ifndef VPTREE_HIERARCHICAL_H_
#define VPTREE_HIERARCHICAL_H_
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
#include <unordered_set>
#include <map>
#include <string>

#include "vptree.cpp"

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "boost/serialization/unordered_map.hpp"
#include "metric_trees_helpers.h"

using namespace Rcpp;
using namespace std;
using namespace boost;

template<typename T>
   class VpTreeHierarchical : public VpTree<T>
{
   typedef typename VpTreeHierarchical<T>::Node Node;
   typedef typename VpTree<T>::HeapItem HeapItem;
   public:
   static const char* ClassName;
   VpTreeHierarchical(int m=4, int minM=2, size_t max_leaf_size=25, size_t vantage_point_candidates=5, size_t test_point_count=15)
      : VpTree<T>(m, minM,  max_leaf_size, vantage_point_candidates, test_point_count)
   {
   }

   ~VpTreeHierarchical() {
   }


   NumericMatrix hierarchicalClustering()
   {
      NumericMatrix ret(this->_items.size()-1, 2);
      neighborsCount = vector<int>(this->_items.size(), 0);
      nearestNeighbors = vector<queue<HeapItem>>(this->_items.size());

      priority_queue<HeapHierarchicalItem> pq;
      //Rcout << "_items.size() = " << this->_items.size() << endl;
      for(int i=0;i<this->_items.size();i++)
      {
         //Rcout << "kazdemu na poczatek znajduje sasiada..." << endl;
         //stop("kazdemu na poczatek znajduje sasiada");
         HeapItem* hi;
         getNearestNeighbor(i, hi);

         if(hi != NULL)
         {
            Rcout <<"dla " << i << "najblizszym jest " << hi->index << endl;
            pq.push(HeapHierarchicalItem(i, hi->index, hi->dist));
         }
         //stop("po pierwszym wstepnym");
      }
      //stop("po wstepnym zebraniu sasiadow");
      std::map<int,int> rank;
      std::map<int,int> parent;

      boost::disjoint_sets<
         associative_property_map<std::map<int,int>>,
         associative_property_map<std::map<int,int>> > ds(
            make_assoc_property_map(rank),
            make_assoc_property_map(parent));

      for(int i=0; i<this->_items.size(); i++)
          ds.make_set(i);
      //stop("po stworzeniu union find");
      //Rcout << "pq size = " << pq.size() << endl;
      //stop(std::to_string(pq.size()));
      int i = 0;
      while(i < this->_items.size() - 1)
      //for(int i=0;i<this->_items.size() - 1 ; i++)
      {
         HeapHierarchicalItem hhi = pq.top();
         pq.pop();

         int s1 = ds.find_set(hhi.index1);
         int s2 = ds.find_set(hhi.index2);
         if(s1 != s2)
         {
            ret(i,0)=hhi.index1;
            ret(i,1)=hhi.index2;
            Rcout << "el1="<<ret(i,0)<< "el2=" <<ret(i,1)<< ", i =" << i << endl;
            i++;
            ds.link(s1, s2);
            Rcout << "el1="<<hhi.index1+1<< "el2=" <<hhi.index2 +1<< endl;

         }
         //stop("przed dociaganiem sasiadow");
         HeapItem* hi;
         getNearestNeighbor(hhi.index1, hi);
         if(hi != NULL)
            pq.push(HeapHierarchicalItem(hhi.index1, hi->index, hi->dist));

         getNearestNeighbor(hhi.index2, hi);
         if(hi != NULL)
            pq.push(HeapHierarchicalItem(hhi.index2, hi->index, hi->dist));
         //stop("po pierwszej iteracji");
      }
      return ret;
   }
   private:
   std::vector<int> neighborsCount;
   std::vector<queue<HeapItem>> nearestNeighbors;

   IntegerMatrix hclust_merge_matrix(const IntegerMatrix x) const {
      // x has 0-based indices
      int n = x.nrow();
      if (x.ncol() != 2) stop("x should have 2 columns");

      IntegerMatrix y(n, 2);
      std::vector< std::unordered_set<int> > curclust(n);

      for (int k=0; k<n; ++k) {
         int i = x(k,0)+1;
         int j = x(k,1)+1;
         int si=k-1, sj=k-1;
         while (si >= 0 && curclust[si].find(i) == curclust[si].end()) si--;
         while (sj >= 0 && curclust[sj].find(j) == curclust[sj].end()) sj--;
         if (si < 0 && sj < 0) {
            curclust[k].insert(i);
            curclust[k].insert(j);
            y(k,0) = -i;
            y(k,1) = -j;
         }
         else if (si < 0 && sj >= 0) {
            curclust[k].insert(curclust[sj].begin(), curclust[sj].end());
            curclust[k].insert(i);
            curclust[sj].clear(); // no longer needed
            y(k,0) = -i;
            y(k,1) = sj+1;
         }
         else if (si >= 0 && sj < 0) {
            curclust[k].insert(curclust[si].begin(), curclust[si].end());
            curclust[k].insert(j);
            curclust[si].clear(); // no longer needed
            y(k,0) = si+1;
            y(k,1) = -j;
         }
         else { // if (si >= 0 && sj >= 0)
            curclust[k].insert(curclust[si].begin(), curclust[si].end());
            curclust[k].insert(curclust[sj].begin(), curclust[sj].end());
            curclust[si].clear(); // no longer needed
            curclust[sj].clear(); // no longer needed
            y(k,0) = si+1;
            y(k,1) = sj+1;
         }
      }
      return y;
   }

   struct HeapHierarchicalItem {
      HeapHierarchicalItem( int index1, int index2, double dist) :
      index1(index1), index2(index2), dist(dist) {}
      int index1;
      int index2;
      double dist;
      bool operator<( const HeapHierarchicalItem& o ) const {
         return dist > o.dist;
      }
   };

   void getNearestNeighbor(int index, HeapItem*& hiret)
   {
      const int delta = 2;//this->_items.size();
      //Rcout << "dociagam sasiadow..." << endl;
      if(nearestNeighbors[index].empty())
      {
         //Rcout << "kolejka pusta, trzeba dociagnac sasiadow..." << endl;
         std::vector<int> results;
         std::vector<double> distances;
         priority_queue<HeapItem> heap;
         searchKNNKnownIndexHierarchical(index, neighborsCount[index]+delta, &results, &distances, false);
         Rcout << "index=" << index << endl;
         for(int i=neighborsCount[index]; i < min(neighborsCount[index]+delta, (int)results.size()); i++)
         {
            Rcout << "results[i]" << results[i] << endl;
            Rcout << "distances[i]" << distances[i] << endl;
            nearestNeighbors[index].push(HeapItem(results[i], distances[i] ));
         }
         //Rcout << "dociagnalem " << heap.size() << "sasiadow" << endl;
         //nearestNeighbors[index] = heap;
         neighborsCount[index] = min(neighborsCount[index]+delta, (int)this->_items.size());
      }

      if(!nearestNeighbors[index].empty())
      {
         hiret = new HeapItem(-1,-1);
         *hiret = nearestNeighbors[index].front();
         nearestNeighbors[index].pop();
      }
      else
      {
         hiret = NULL;
         //stop("nie ma sasiadow!");
      }
   }

   void searchKNNKnownIndexHierarchical(int index, int k, std::vector<int>* results,
                            std::vector<double>* distances, bool findItself = true)
   {
      if(index < 0 || index >= this->_items.size()) stop("Index out of bounds.");
      #ifdef DEBUG
      _distance.metricCalculated = 0;
      _distance.hashmapHit = 0;
      #endif
      std::priority_queue<typename VpTreeHierarchical<T>::HeapItem> heap;

      this->_tau = std::numeric_limits<double>::max();
      search( VpTree<T>::_root, index, true, k, heap, findItself );

      results->clear(); distances->clear();
      Rcout << "znalazlem " << heap.size() << "sasiadow" << endl;
      while( !heap.empty() ) {
         results->push_back( heap.top().index );
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



   void search( Node* node, int index, bool isKNN, int k,
               std::priority_queue<HeapItem>& heap, bool findItself )
   {
      if ( node == NULL ) return;

      //printf("dist=%g tau=%gn", dist, _tau );
      if(node->isLeaf)
      {
         for(size_t i=0;i<node->points.size();i++)
         {
            double dist2 = this->_distance(node->points[i], index );
            if(dist2 < 1e-6)
            {
                if(!findItself && R_compute_identical(this->_items[node->points[i]], this->_items[index], 16))
                   continue;
            }
            if ( (dist2 < this->_tau && isKNN) || (dist2 <= this->_tau && !isKNN) ) {
               if(node->points[i] > index)
               {
                  if ( heap.size() ==(size_t) k && isKNN) heap.pop();
                  heap.push(HeapItem(node->points[i], dist2) );
                  if ( heap.size() == (size_t) k && isKNN)
                  {
                     this->_tau = heap.top().dist;
                     //Rcout << "current tau=" << _tau << endl;
                  }
               }
            }
         }
      }
      else
      {
         double dist = this->_distance(node->vpindex, index);
         vector<bool> visited(node->childCount, false);

         for(int i=0;i<node->childCount-1;i++)
         {
            //if ( dist < node->radiuses[i] ) {
            if ( dist - this->_tau <= node->radiuses[i] && !visited[i]) {
               search( node->children[i], index, isKNN, k, heap, findItself );
               visited[i]=true;
            }
            if ( dist + this->_tau >= node->radiuses[i] && !visited[i+1]) {
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

template<typename T>
   const char* VpTreeHierarchical<T>::ClassName = "VpTreeHierarchical";

template<typename T>
   void checkIsVpTreeHierarchicalClass(XPtr< VpTreeHierarchical<T> >& _tree)
{
   if(strcmp(_tree.attr("class"), VpTreeHierarchical<T>::ClassName))
      stop("not a VpTreeHierarchical object");
}

// [[Rcpp::export]]
NumericMatrix hierarchical(Function distance, List listobj,
                   bool isSimilarity = false,
                   int m=2,
                   int minm=4,
                   int maxLeafPointsCount = 25,
                   int vantageCandidatesCount = 5,
                   int testPointsCount = 15) { //https://code.google.com/p/vptree/source/browse/src/vptree/VpTreeNode.java
   VpTreeHierarchical<RObject>* _tree = new VpTreeHierarchical<RObject>(m, minm, maxLeafPointsCount, vantageCandidatesCount, testPointsCount);
   distClass distci;
   distci.isSimilarity = isSimilarity;
   distci.distance = new RFunction(distance);
   _tree->setDistanceFunction(distci);

   vector<RObject> points = createStdVectorOfRobjects(listobj);
   //stop("przed kreacja");
   Rcout << "przed kreacja2" << endl;
   _tree->create(points);
   //stop("przed klasteringiem");
   return _tree->hierarchicalClustering();
}
#endif /* VPTREE_HIERARCHICAL_H_ */
