#ifndef __HCLUST2_DISTANCE_H
#define __HCLUST2_DISTANCE_H

// ---------------------------------------------------------------------
#define HASHMAP_DISABLE
#define HASHMAP_COUNTERS
#define VERBOSE 0
// ---------------------------------------------------------------------



#include <boost/functional/hash.hpp>
#include <vector>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <Rcpp.h>



namespace DataStructures{

struct SortedPoint
{
   size_t i;
   size_t j;

   SortedPoint()
      :i(0),j(0) {}

   SortedPoint(size_t _i, size_t _j)
   {
      if(_j < _i)
      {
         i = _j;
         j = _i;
      }
      else
      {
         i = _i;
         j = _j;
      }
   }

   bool operator==(const SortedPoint &other) const
   {
      return (i == other.i && j == other.j);
   }
};

} // namespace DataStructures


namespace std {

   template <>
      struct hash<DataStructures::SortedPoint>
   {
      std::size_t operator()(const DataStructures::SortedPoint& k) const
      {
        std::size_t seed = 0;
        boost::hash_combine(seed, k.i);
        boost::hash_combine(seed, k.j);
        return seed;
      }
   };
} // namespace std


namespace DataStructures {

struct Distance {
#ifndef HASHMAP_DISABLE
   std::vector< std::unordered_map<size_t, double> > hashmap;
#ifdef HASHMAP_COUNTERS
   size_t hashmapHit;
   size_t hashmapMiss;
#endif
#endif
   size_t n;

   Distance(size_t n);
   virtual ~Distance();
   inline size_t getObjectCount() { return n; }
   virtual double compute(size_t v1, size_t v2)  const = 0;
   static Distance* createDistance(Rcpp::RObject objects, Rcpp::RObject distance);

#ifndef HASHMAP_DISABLE
   double operator()(size_t v1, size_t v2);
#else
   inline double operator()(size_t v1, size_t v2)  const {
      return compute(v1, v2);
   }
#endif
};

struct EuclideanDistance : public Distance
{
   SEXP robj1;
   const double* items;
   size_t m;

   EuclideanDistance(const Rcpp::NumericMatrix& points) :
      Distance(points.nrow()), robj1(points),
      items(REAL((SEXP)points)), m(points.ncol())
   {
      R_PreserveObject(robj1);
   }

#ifdef HASHMAP_COUNTERS
   virtual ~EuclideanDistance()
   {
#if VERBOSE > 5
      Rprintf("[%010.3f] destroying distance object\n", clock()/(float)CLOCKS_PER_SEC);
#endif
      R_ReleaseObject(robj1);
   }
#endif

   virtual double compute(size_t v1, size_t v2) const;
};


struct ManhattanDistance : public Distance
{
   SEXP robj1;
   const double* items;
   size_t m;

   ManhattanDistance(const Rcpp::NumericMatrix& points) :
      Distance(points.nrow()), robj1(points),
      items(REAL((SEXP)points)), m(points.ncol())
   {
      R_PreserveObject(robj1);
   }

   virtual ~ManhattanDistance()
   {
#if VERBOSE > 5
      Rprintf("[%010.3f] destroying distance object\n", clock()/(float)CLOCKS_PER_SEC);
#endif
      R_ReleaseObject(robj1);
   }

   double compute(size_t v1, size_t v2) const;
};


struct MaximumDistance : public Distance
{
   SEXP robj1;
   const double* items;
   size_t m;
#ifdef HASHMAP_COUNTERS
   size_t hashmapHit;
   size_t hashmapMiss;
#endif

   MaximumDistance(const Rcpp::NumericMatrix& points) :
      Distance(points.nrow()), robj1(points),
      items(REAL((SEXP)points)), m(points.ncol())
   {
      R_PreserveObject(robj1);
   }

   virtual ~MaximumDistance()
   {
#if VERBOSE > 5
      Rprintf("[%010.3f] destroying distance object\n", clock()/(float)CLOCKS_PER_SEC);
#endif
      R_ReleaseObject(robj1);
   }

   virtual double compute(size_t v1, size_t v2) const;
};

struct GenericRDistance : public Distance
{
   Rcpp::Function distfun;
   std::vector<Rcpp::RObject> items;

   GenericRDistance(const Rcpp::Function& distfun, const std::vector<Rcpp::RObject>& items) :
      Distance(items.size()),
      distfun(distfun),
      items(items)
   {
      R_PreserveObject(distfun);
   }

   virtual ~GenericRDistance()
   {
      R_ReleaseObject(distfun);
   }

   virtual double compute(size_t v1, size_t v2) const;
};

};




#endif
