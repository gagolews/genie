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
#include "mtree/mtree.h"

using namespace Rcpp;
using namespace std;
using namespace mt;

template<typename Data, typename DistanceFunction, typename SplitFunction>
   const char* mtree<Data, DistanceFunction, SplitFunction>::ClassName = "mTree";

template<typename Data, typename DistanceFunction, typename SplitFunction>
   void checkIsMTreeClass(XPtr< mtree<Data, DistanceFunction, SplitFunction> >& _tree)
{
   if(strcmp(_tree.attr("class"), mtree<Data, DistanceFunction, SplitFunction>::ClassName))
      stop("not an mTree object");
}


//' @rdname mtree
//' @details
//' \code{mtree_create} creates an empty m-tree instance. Please use
//' \code{mtree_build} for populating the tree.
//'
//' @return
//' \code{mtree_create} returns a new, empty m-tree.
//' @param distance an R function, which calculates metric in the space
//' @param isSimilarity a logical value, does a distance function calculate
//' distance or similarity?
//' @param m single integer, maximum number of children can be in a node
//' in the tree
//' @param minm single integer, minimum number of children can be in a node
//' in the tree
//' @param maxLeafPointsCount single integer, maximum number of
//' points in a single leaf
//' @param vantageCandidatesCount single integer, number of candidates
//' considered as a vantage point in a single node. Too big number can
//' cause performance problem in building the tree, but too small can cause
//' performance problem when searching.
//' @param testPointsCount single integer, how many points are considered with
//' a single candidate for vantage point in a node to assess a variance of
//' distances
// [[Rcpp::export]]
SEXP mtree_create(Function distance, bool isSimilarity = false) {
   distClass distci;
   distci.isSimilarity = isSimilarity;
   distci.distance = new RFunction(distance);
   mtree<RObject,distClass>* vec = new mtree<RObject,distClass>(mtree<RObject,distClass>::DEFAULT_MIN_NODE_CAPACITY, -1, distci);
   XPtr< mtree<RObject,distClass> > retval =  XPtr< mtree<RObject,distClass> >(vec, true);
   retval.attr("class") = mtree<RObject,distClass>::ClassName;
   return retval;
}

//' @rdname mtree
//' @details
//' \code{mtree_insert} inserts a new point into a m-tree.
//' Calling this function n times has worse performance in comparison
//' to calling one \code{mtree_build} with n points in the beginning.
//'
//' @return
//' \code{mtree_insert} does not return anything interesting.
//' @param obj a point to insert
// [[Rcpp::export]]
void mtree_insert(SEXP tree, RObject obj) {
   XPtr< mtree<RObject,distClass> > _tree = Rcpp::as< XPtr< mtree<RObject,distClass> > > (tree);
   checkIsMTreeClass(_tree);
   R_PreserveObject(obj);
   (*_tree).add(obj);
}



//' @rdname mtree
//' @details
//' \code{mtree_build} populates an empty tree with set of points. If called
//' for non-empty tree, whole structure of previous tree is deleted. This is
//' recommended way of building a tree (instead of inserts by
//' \code{mtree_insert}).
//'
//' @return
//' \code{mtree_build} does not return anything interesting.
//' @param listobj an R list, where every element of list is assumed as one
//' point in a  metric space
// [[Rcpp::export]]
void mtree_build(SEXP tree, List listobj) {
   XPtr< mtree<RObject,distClass> > _tree = Rcpp::as< XPtr< mtree<RObject,distClass> > > (tree);
   checkIsMTreeClass(_tree);
   vector<RObject> points = createStdVectorOfRobjects(listobj);
   for(int i=0;i<points.size();i++)
      (*_tree).add(points[i]);
}

//' @rdname mtree
//' @details
//' \code{mtree_searchKNN} finds k nearest neighbours for a given object
//' in the tree.
//' Please use this function for objects, which are not in a set of points
//' inserted to the tree or which were used for building the tree. For these
//' cases please use \code{mtree_searchKNNKnown} or
//' \code{mtree_searchKNNKnownIndex} which can be faster.
//'
//' @return
//' \code{mtree_searchKNN} returns a list. First element of the list
//' is a list with found elements. A second element of the list is a list
//' with distances these objects from a given object.
//' @param p an R object for which neighbours are found
//' @param k a single integer, number of neighbours to find
//' @param findItself boolean value, should results contain an given object?
// [[Rcpp::export]]
List mtree_searchKNN(SEXP tree, RObject p, int k, bool findItself = true)
{
   XPtr< mtree<RObject,distClass> > _tree = Rcpp::as< XPtr< mtree<RObject,distClass> > > (tree);
   checkIsMTreeClass(_tree);
   auto query = (*_tree).get_nearest_by_limit(p, k);
   vector<mtree<RObject,distClass>::query::result_item> results(query.begin(), query.end());

   size_t resultsSize = results.size();
   if(!findItself)
   {
      for(int i=0;i<results.size();i++)
      {
          if(R_compute_identical(results[i].data, p, 16))
          {
              resultsSize--;
          }
      }
   }

   List resultList(2);
   List elements(resultsSize);
   List distancesRcpp(resultsSize);
   int index = 0;
   for(int i = 0; i < results.size(); ++i)
   {
      if(!findItself && R_compute_identical(results[i].data, p, 16)) continue;
      elements[index] = results[i].data;
      distancesRcpp[index] = results[i].distance;
      index++;
   }
   resultList[0] = elements;
   resultList[1] = distancesRcpp;
   return resultList;
}

//' @rdname mtree
//' @details
//' \code{mtree_searchRadius} finds neighbours for a given object
//' within a given radius in the tree.
//' Please use this function for objects, which are not in a set of points
//' inserted to the tree or which were used for building the tree. For these
//' cases please use \code{mtree_searchRadiusKnown} or
//' \code{mtree_searchRadiusKnownIndex} which can be faster.
//'
//' @return
//' \code{mtree_searchRadius} returns a list. First element of the list
//' is a list with found elements. A second element of the list is a list
//' with distances these objects from a given object.
//' @param p an R object for which neighbours are found
//' @param tau a float value, a radius
//' @param findItself boolean value, should results contain an given object?
// [[Rcpp::export]]
List mtree_searchRadius(SEXP tree, RObject p, double tau, bool findItself = true)
{
   XPtr< mtree<RObject,distClass> > _tree = Rcpp::as< XPtr< mtree<RObject,distClass> > > (tree);
   checkIsMTreeClass(_tree);
   auto query = (*_tree).get_nearest(p, tau, std::numeric_limits<unsigned int>::max());
   vector<mtree<RObject,distClass>::query::result_item> results(query.begin(), query.end());

  size_t resultsSize = results.size();
   if(!findItself)
   {
      for(int i=0;i<results.size();i++)
      {
          if(R_compute_identical(results[i].data, p, 16))
          {
              resultsSize--;
          }
      }
   }

   List resultList(2);
   List elements(resultsSize);
   List distancesRcpp(resultsSize);
   int index = 0;
   for(int i = 0; i < results.size(); ++i)
   {
      if(!findItself && R_compute_identical(results[i].data, p, 16)) continue;
      elements[index] = results[i].data;
      distancesRcpp[index] = results[i].distance;
      index++;
   }
   resultList[0] = elements;
   resultList[1] = distancesRcpp;
   return resultList;
}
