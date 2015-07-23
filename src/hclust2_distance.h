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

class Distance {
private:
#ifndef HASHMAP_DISABLE
   std::vector< std::unordered_map<size_t, double> > hashmap;
#ifdef HASHMAP_COUNTERS
   size_t hashmapHit;
   size_t hashmapMiss;
#endif
#endif
protected:
   size_t n;
   virtual double compute(size_t v1, size_t v2)  const = 0;

public:
   Distance(size_t n);
   virtual ~Distance();
   inline size_t getObjectCount() { return n; }
   static Distance* createDistance(Rcpp::RObject objects, Rcpp::RObject distance);

#ifndef HASHMAP_DISABLE
   double operator()(size_t v1, size_t v2);
#else
   inline double operator()(size_t v1, size_t v2)  const {
      return compute(v1, v2);
   }
#endif
};


class GenericMatrixDistance : public Distance
{
protected:
   double* items;
   size_t m;

public:
   GenericMatrixDistance(const Rcpp::NumericMatrix& points) :
      Distance(points.nrow()),
      items(REAL((SEXP)points)), m(points.ncol())
   {
      // act on a transposed matrix to avoid many L1/L... cache misses
      items = new double[m*n];
      const double* items2 = REAL((SEXP)points);
      double* items_ptr = items;
      for (size_t i=0; i<n; ++i)
         for (size_t j=0; j<m; ++j)
            *(items_ptr++) = items2[j*n+i];
   }

   virtual ~GenericMatrixDistance()
   {
#if VERBOSE > 5
      Rprintf("[%010.3f] destroying distance object\n", clock()/(float)CLOCKS_PER_SEC);
#endif
      delete [] items;
   }
};


class EuclideanDistance : public GenericMatrixDistance
{
// private:
//    std::vector<double> sqobs;

protected:
   virtual double compute(size_t v1, size_t v2) const;

public:
   EuclideanDistance(const Rcpp::NumericMatrix& points) :
      GenericMatrixDistance(points)
   {
//       const double* items_ptr = items;
//       for (size_t i=0; i<n; ++i) {
//          double sqobs_cur = 0.0;
//          for (size_t j=0; j<m; ++j) {
//             sqobs_cur += (*items_ptr)*(*items_ptr);
//             ++items_ptr;
//          }
//          sqobs[i] = sqobs_cur*0.5;
//       }
   }
};


class ManhattanDistance : public GenericMatrixDistance
{
protected:
   double compute(size_t v1, size_t v2) const;

public:
   ManhattanDistance(const Rcpp::NumericMatrix& points) :
      GenericMatrixDistance(points)  {   }
};


class MaximumDistance : public GenericMatrixDistance
{
protected:
   double compute(size_t v1, size_t v2) const;

public:
   MaximumDistance(const Rcpp::NumericMatrix& points) :
      GenericMatrixDistance(points)  {   }
};


class GenericRDistance : public Distance
{
private:
   Rcpp::Function distfun;
   std::vector<Rcpp::RObject> items;

protected:
   virtual double compute(size_t v1, size_t v2) const;

public:
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
};


class DistObjectDistance : public Distance
{
protected:
   SEXP robj1;
   const double* items;

protected:
   virtual double compute(size_t v1, size_t v2) const;

public:
   DistObjectDistance(const Rcpp::NumericVector& distobj) :
      Distance((size_t)((Rcpp::NumericVector)distobj.attr("Size"))[0]),
      robj1(distobj),
      items(REAL((SEXP)distobj))
   {
      if ((size_t)XLENGTH((SEXP)distobj) != n*(n-1)/2)
         Rcpp::stop("incorrect dist object length.");
      R_PreserveObject(robj1);
   }

   virtual ~DistObjectDistance()
   {
      R_ReleaseObject(robj1);
   }
};

} // namespace DataStructures

#endif
