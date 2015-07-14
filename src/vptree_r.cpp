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

#include "vptree.cpp"

//' @rdname vptree
//' @details
//' \code{vptree_create} creates an empty vp-tree instance. Please use
//' \code{vptree_build} for populating the tree.
//'
//' @return
//' \code{vptree_create} returns a new, empty vp-tree.
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
SEXP vptree_create(Function distance,
                   bool isSimilarity = false,
                   int m=2,
                   int minm=4,
                   int maxLeafPointsCount = 25,
                   int vantageCandidatesCount = 5,
                   int testPointsCount = 15) { //https://code.google.com/p/vptree/source/browse/src/vptree/VpTreeNode.java
   VpTree<RObject>* vec = new VpTree<RObject>(m, minm, maxLeafPointsCount, vantageCandidatesCount, testPointsCount);
   XPtr< VpTree<RObject> > retval =  XPtr< VpTree<RObject> >(vec, true);
   retval.attr("class") = VpTree<RObject>::ClassName;
   distClass distci;
   distci.isSimilarity = isSimilarity;
   distci.distance = new RFunction(distance);
   (*retval).setDistanceFunction(distci);

   return retval;
}

//' @rdname vptree
//' @details
//' \code{vptree_insert} inserts a new point into a vp-tree.
//' Calling this function n times has worse performance in comparison
//' to calling one \code{vptree_build} with n points in the beginning.
//'
//' @return
//' \code{vptree_insert} does not return anything interesting.
//' @param obj a point to insert
// [[Rcpp::export]]
void vptree_insert(SEXP tree, RObject obj) {
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   R_PreserveObject(obj);
   (*_tree).insert(obj);
}

//' @rdname vptree
//' @details
//' \code{vptree_set_distancefunction} sets a metric function for a tree.
//' Setting this function after building a tree will cause problems. Do not
//' use if you do not know what you are doing.
//'
//' @return
//' \code{vptree_set_distancefunction} does not return anything interesting.
//' @param distance an R function, which calculates metric in the space
//' @param isSimilarity a logical value, does a distance function calculate
//' distance or similarity?
// [[Rcpp::export]]
void vptree_set_distancefunction(SEXP tree, Function distance, bool isSimilarity = false) {
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   distClass distci;
   distci.isSimilarity = isSimilarity;
   distci.distance = new RFunction(distance);
   (*_tree).setDistanceFunction(distci);
}

//' @rdname vptree
//' @details
//' \code{vptree_build} populates an empty tree with set of points. If called
//' for non-empty tree, whole structure of previous tree is deleted. This is
//' recommended way of building a tree (instead of inserts by
//' \code{vptree_insert}).
//'
//' @return
//' \code{vptree_build} does not return anything interesting.
//' @param listobj an R list, where every element of list is assumed as one
//' point in a  metric space
// [[Rcpp::export]]
void vptree_build(SEXP tree, List listobj) {
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   vector<RObject> points = createStdVectorOfRobjects(listobj);
   (*_tree).create(points);
}

//' @rdname vptree
//' @details
//' \code{vptree_searchKNN} finds k nearest neighbours for a given object
//' in the tree.
//' Please use this function for objects, which are not in a set of points
//' inserted to the tree or which were used for building the tree. For these
//' cases please use \code{vptree_searchKNNKnown} or
//' \code{vptree_searchKNNKnownIndex} which can be faster.
//'
//' @return
//' \code{vptree_searchKNN} returns a list. First element of the list
//' is a list with found elements. A second element of the list is a list
//' with distances these objects from a given object.
//' @param p an R object for which neighbours are found
//' @param k a single integer, number of neighbours to find
//' @param findItself boolean value, should results contain an given object?
// [[Rcpp::export]]
List vptree_searchKNN(SEXP tree, RObject p, int k, bool findItself = true)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   std::vector<RObject> results;
   std::vector<double> distances;
   (*_tree).searchKNN( p, k, &results, &distances, findItself);

   List resultList(2);
   List elements(results.size());
   List distancesRcpp(results.size());
   int index =0;
   for (std::vector<RObject>::iterator it = results.begin(); it!=results.end(); ++it)
   {
      elements[index] = *it;
      distancesRcpp[index] = distances[index];
      index++;
   }
   resultList[0] = elements;
   resultList[1] = distancesRcpp;
   return resultList;
}

//' @rdname vptree
//' @details
//' \code{vptree_searchKNNKnown} finds k nearest neighbours for a given object
//' in the tree which was used for building the tree or inserted after that.
//' Please do not use this function for objects, which are "new" in a metric
//' space. This function can be faster than \code{vptree_searchKNN}.
//'
//' @return
//' \code{vptree_searchKNNKnown} returns a list. First element of the list
//' is a list with found elements. A second element of the list is a list
//' with distances these objects from a given object.
//' @param p an R object for which neighbours are found
//' @param k a single integer, number of neighbours to find
//' @param findItself boolean value, should results contain an given object?
// [[Rcpp::export]]
List vptree_searchKNNKnown(SEXP tree, RObject p, int k, bool findItself = true)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   std::vector<RObject> results;
   std::vector<double> distances;
   (*_tree).searchKNNKnown( p, k, &results, &distances, findItself);

   List resultList(2);
   List elements(results.size());
   List distancesRcpp(results.size());
   int index =0;
   for (std::vector<RObject>::iterator it = results.begin(); it!=results.end(); ++it)
   {
      elements[index] = *it;
      distancesRcpp[index] = distances[index];
      index++;
   }
   resultList[0] = elements;
   resultList[1] = distancesRcpp;
   return resultList;
}

//' @rdname vptree
//' @details
//' \code{vptree_searchKNN} finds k nearest neighbours for a given index of
//' the object in space. Indices are as in a list given for building a tree.
//' For inserted objects they are pushed back at the end.
//' Please do not use this function for objects, which are "new" in a metric
//' space. This function can be faster than \code{vptree_searchKNN} and
//' \code{vptree_searchKNNKnown}. Indexing starts with 1, as in R language.
//'
//' @return
//' \code{vptree_searchKNNKnown} returns a list. First element of the list
//' is a list with found elements. A second element of the list is a list
//' with distances these objects from a given object.
//' @param index index of an object in the space
//' @param k a single integer, number of neighbours to find
//' @param findItself boolean value, should results contain an given object?
// [[Rcpp::export]]
List vptree_searchKNNKnownIndex(SEXP tree, int index, int k, bool findItself = true)
{
   index--;
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   std::vector<RObject> results;
   std::vector<double> distances;
   (*_tree).searchKNNKnownIndex( index, k, &results, &distances, findItself);

   List resultList(2);
   List elements(results.size());
   List distancesRcpp(results.size());
   int idx =0;
   for (std::vector<RObject>::iterator it = results.begin(); it!=results.end(); ++it)
   {
      elements[idx] = *it;
      distancesRcpp[idx] = distances[idx];
      idx++;
   }
   resultList[0] = elements;
   resultList[1] = distancesRcpp;
   return resultList;
}

//' @rdname vptree
//' @details
//' \code{vptree_searchRadius} finds neighbours for a given object
//' within a given radius in the tree.
//' Please use this function for objects, which are not in a set of points
//' inserted to the tree or which were used for building the tree. For these
//' cases please use \code{vptree_searchRadiusKnown} or
//' \code{vptree_searchRadiusKnownIndex} which can be faster.
//'
//' @return
//' \code{vptree_searchRadius} returns a list. First element of the list
//' is a list with found elements. A second element of the list is a list
//' with distances these objects from a given object.
//' @param p an R object for which neighbours are found
//' @param tau a float value, a radius
//' @param findItself boolean value, should results contain an given object?
// [[Rcpp::export]]
List vptree_searchRadius(SEXP tree, RObject p, double tau, bool findItself = true)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   std::vector<RObject> results;
   std::vector<double> distances;
   (*_tree).searchRadius( p, tau, &results, &distances, findItself);

   List resultList(2);
   List elements(results.size());
   List distancesRcpp(results.size());
   int idx =0;
   for (std::vector<RObject>::iterator it = results.begin(); it!=results.end(); ++it)
   {
      elements[idx] = *it;
      distancesRcpp[idx] = distances[idx];
      idx++;
   }
   resultList[0] = elements;
   resultList[1] = distancesRcpp;
   return resultList;
}

//' @rdname vptree
//' @details
//' \code{vptree_searchRadiusKnown} finds nearest neighbours for a given object
//' within a given radius
//' in the tree which was used for building the tree or inserted after that.
//' Please do not use this function for objects, which are "new" in a metric
//' space. This function can be faster than \code{vptree_searchRadius}.
//'
//' @return
//' \code{vptree_searchRadiusKnown} returns a list. First element of the list
//' is a list with found elements. A second element of the list is a list
//' with distances these objects from a given object.
//' @param p an R object for which neighbours are found
//' @param tau a float value, a radius
//' @param findItself boolean value, should results contain an given object?
// [[Rcpp::export]]
List vptree_searchRadiusKnown(SEXP tree, RObject p, double tau, bool findItself = true)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   std::vector<RObject> results;
   std::vector<double> distances;
   (*_tree).searchRadiusKnown( p, tau, &results, &distances, findItself);

   List resultList(2);
   List elements(results.size());
   List distancesRcpp(results.size());
   int idx =0;
   for (std::vector<RObject>::iterator it = results.begin(); it!=results.end(); ++it)
   {
      elements[idx] = *it;
      distancesRcpp[idx] = distances[idx];
      idx++;
   }
   resultList[0] = elements;
   resultList[1] = distancesRcpp;
   return resultList;
}

//' @rdname vptree
//' @details
//' \code{vptree_searchRadiusKnownIndex} finds nearest neighbours
//' for a given index of the object in space
//' within a given radius . Indices are as in a list given for building a tree.
//' For inserted objects they are pushed back at the end.
//' Please do not use this function for objects, which are "new" in a metric
//' space. This function can be faster than \code{vptree_searchRadius} and
//' \code{vptree_searchRadiusKnown}. Indexing starts with 1, as in R language.
//'
//' @return
//' \code{vptree_searchRadiusKnownIndex} returns a list. First element of the list
//' is a list with found elements. A second element of the list is a list
//' with distances these objects from a given object.
//' @param index index of an object in the space
//' @param tau a float value, a radius
//' @param findItself boolean value, should results contain an given object?
// [[Rcpp::export]]
List vptree_searchRadiusKnownIndex(SEXP tree, int index, double tau, bool findItself = true)
{
   index--;
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   std::vector<RObject> results;
   std::vector<double> distances;
   (*_tree).searchRadiusKnownIndex( index, tau, &results, &distances, findItself);

   List resultList(2);
   List elements(results.size());
   List distancesRcpp(results.size());
   int idx =0;
   for (std::vector<RObject>::iterator it = results.begin(); it!=results.end(); ++it)
   {
      elements[idx] = *it;
      distancesRcpp[idx] = distances[idx];
      idx++;
   }
   resultList[0] = elements;
   resultList[1] = distancesRcpp;
   return resultList;
}

//' @rdname vptree
//' @details
//' \code{vptree_getItems} returns a list of points in it.
//'
//' @return
//' \code{vptree_getItems} returns a list.
// [[Rcpp::export]]
List vptree_getItems(SEXP vptree)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (vptree);
   checkIsVpTreeClass(_tree);
   return wrap(_tree->getCopyItems());
   //return List();
}

//' @rdname vptree
//' @details
//' \code{vptree_getFunction} returns a metric function used in the tree.
//'
//' @return
//' \code{vptree_getFunction} returns a function.
// [[Rcpp::export]]
Function vptree_getFunction(SEXP vptree)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (vptree);
   checkIsVpTreeClass(_tree);
   return _tree->getMetricFunction();
}

//' @rdname vptree
//' @details
//' \code{vptree_serialize} saves a structure of a tree to a file.
//' This function do not save a list of points or metric function. For saving whole
//' tree please use \code{vptree_save}. Do not use this function unless you know what
//' you are doing.
//'
//' @return
//' \code{vptree_getFunction} does not return anything interesting.
//' @param filename a path to file in which a strucutre of a tree will be saved.
// [[Rcpp::export]]
void vptree_serialize(SEXP vptree, std::string filename)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (vptree);
   checkIsVpTreeClass(_tree);

   std::ofstream ofs(filename);
   {
      if(!ofs.good())
      {
         stop("File creation error.");
      }
      // save tree to archive
      boost::archive::text_oarchive oa(ofs);
      // write class instance to archive
      oa << *_tree;
   }
}

//' @rdname vptree
//' @details
//' \code{vptree_read} loads a structure of a tree from a file.
//' This function do not load a list of points or metric function. For loading whole
//' tree please use \code{vptree_load}. Do not use this function unless you know what
//' you are doing.
//'
//' @return
//' \code{vptree_read} does not return anything interesting.
//' @param filename a path to file in which a strucutre of a tree is saved.
// [[Rcpp::export]]
SEXP vptree_read(std::string filename)
{
   VpTree<RObject>* vptree_read = new VpTree<RObject>(4, 2, 25, 5, 10);
   // create and open an archive for input
   {
      std::ifstream ifs(filename);
      boost::archive::text_iarchive ia(ifs);
      // read class state from archive
      ia >> *vptree_read;
   }
   XPtr< VpTree<RObject> > retval =  XPtr< VpTree<RObject> >(vptree_read, true);
   retval.attr("class") = VpTree<RObject>::ClassName;
   return retval;
}

//' @rdname vptree
//' @details
//' \code{vptree_setItems} sets a set of points in metric space. Using this
//' function after building a tree will cause problems. Do not use this function
//' unless you know what you are doing.
//'
//' @return
//' \code{vptree_setItems} does not return anything interesting.
//' @param items a list of points.
// [[Rcpp::export]]
void vptree_setItems(SEXP vptree, List items)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (vptree);
   checkIsVpTreeClass(_tree);
   vector<RObject> points = createStdVectorOfRobjects(items);
   _tree->setItems(points);
}

//' @rdname vptree
//' @details
//' \code{vptree_setMetricFunction} sets metric function. Using this
//' function after building a tree will cause problems. Do not use this function
//' unless you know what you are doing.
//'
//' @return
//' \code{vptree_setMetricFunction} does not return anything interesting.
//' @param f a metric function.
// [[Rcpp::export]]
void vptree_setMetricFunction(SEXP vptree, Function f)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (vptree);
   checkIsVpTreeClass(_tree);
   _tree->setMetricFunction(f);
}

//' @rdname vptree
//' @details
//' \code{vptree_treeSize} returns the size of a tree. To returned value adding
//' a size of RObjects in space is needed. In other words, this function
//' returns size of all nodes in tree (and a tree itself),
//' including indices of points in space,
//' size of pointers, radiuses etc., but no points (RObjects) underlying.
//'
//'
//' @return
//' \code{vptree_treeSize} returns a size in bytes of a vp-tree.
// [[Rcpp::export]]
size_t vptree_treeSize(SEXP vptree)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (vptree);
   checkIsVpTreeClass(_tree);
   return _tree->treeSize();
}

//' @rdname vptree
//' @details
//' \code{vptree_treeSize} returns the height of a tree.
//'
//' @return
//' \code{vptree_treeSize} returns the height of a tree.
// [[Rcpp::export]]
int vptree_treeHeight(SEXP vptree)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (vptree);
   checkIsVpTreeClass(_tree);
   return _tree->treeHeight();
}

#ifdef DEBUG
//' @rdname vptree
//' @details
//' \code{vptree_printCounters} prints how many hits every element in hashmap got.
//'
//' @return
//' \code{vptree_printCounters} does not return anything interesting.
// [[Rcpp::export]]
void vptree_printCounters(SEXP vptree)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (vptree);
   checkIsVpTreeClass(_tree);
   _tree->printCounters();
}
#else
void vptree_printCounters(SEXP vptree)
{
   stop("Please compile with DEBUG flag on, e.g. -DDEBUG.");
}
#endif
