/* ************************************************************************* *
 *   This file is part of the `DataStructures` package.                      *
 *                                                                           *
 *   Copyright 2015 Maciej Bartoszuk, Anna Cena, Marek Gagolewski,           *
 *                                                                           *
 *   Parts of the code are taken from the 'CITAN' R package by M. Gagolewski *
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


#ifndef __HCLUST2_DISTANCE_H
#define __HCLUST2_DISTANCE_H

// ---------------------------------------------------------------------
// #define HASHMAP_ENABLED
#define GENERATE_STATS
#define VERBOSE 1
// ---------------------------------------------------------------------



/* to do: dist for CharacterVector (objects=strings)
   dists = levensthein (q-gram: not -> see matrix input on q-gram profiles), lcs, dam-lev

  numeric -> metric: hamming, binary (see dist) minkowski (p), canberra

  allow external ptr distance:

 double dist(const char* s1, int nx, const char* s2, int ny)

 double dist(const double* s1, int nx, const double* s2, int ny)

 double dist(const int* s1, int nx, const int* s2, int ny)

 double dist(SEXP s1, SEXP s2)
*/


#include <boost/functional/hash.hpp>
#include <vector>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <Rcpp.h>


namespace DataStructures {

struct DistanceStats {
   size_t hashmapHit;
   size_t hashmapMiss;
   size_t distCallCount;
   size_t distCallTheoretical;

   DistanceStats(size_t n) :
      hashmapHit(0), hashmapMiss(0), distCallCount(0),
      distCallTheoretical(n*(n-1)/2) {}

   ~DistanceStats() {
#if VERBOSE > 0
#if defined(HASHMAP_ENABLED) && defined(GENERATE_STATS)
   Rprintf("             distance function hashmap #hits: %.0f, #miss: %.0f, est.mem.used: ~%.1fMB (vs %.1fMB)\n",
      (double)hashmapHit, (double)hashmapMiss,
      8.0f*hashmapMiss/1000.0f/1000.0f,
      8.0f*distCallTheoretical/1000.0f/1000.0f);
#endif
#if defined(GENERATE_STATS)
   Rprintf("             distance function total calls: %.0f (i.e., %.2f%% of %.0f)\n",
      (double)distCallCount,
      (double)distCallCount*100.0/(double)distCallTheoretical,
      (double)distCallTheoretical
   );
#endif
#endif
   }

   Rcpp::NumericVector toR() const {
      return Rcpp::NumericVector::create(
         Rcpp::_["hashmapHit"]
            = (hashmapHit>0)?(double)hashmapHit:NA_REAL,
         Rcpp::_["hashmapMiss"]
            = (hashmapMiss>0)?(double)hashmapMiss:NA_REAL,
         Rcpp::_["distCallCount"]
            = (distCallCount>0)?(double)distCallCount:NA_REAL,
         Rcpp::_["distCallTheoretical"]
            = (distCallTheoretical>0)?(double)distCallTheoretical:NA_REAL
      );
   }
};


class Distance {
private:
#ifdef HASHMAP_ENABLED
   std::vector< std::unordered_map<size_t, double> > hashmap;
#endif
   DistanceStats stats;

protected:
   size_t n;
   virtual double compute(size_t v1, size_t v2)  const = 0;

public:
   Distance(size_t n);
   virtual ~Distance();
   inline size_t getObjectCount() { return n; }
   static Distance* createDistance(Rcpp::RObject distance, Rcpp::RObject objects);

   inline const DistanceStats& getStats() { return stats; }

#ifdef HASHMAP_ENABLED
   double operator()(size_t v1, size_t v2);
#else
#ifndef GENERATE_STATS
#define GENERATE_STATS_CONST const
#else
#define GENERATE_STATS_CONST /* const */
#endif
   inline double operator()(size_t v1, size_t v2) GENERATE_STATS_CONST {
#ifdef GENERATE_STATS
      ++stats.distCallCount;
#endif
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
// #if VERBOSE > 5
//       Rprintf("[%010.3f] destroying distance object\n", clock()/(float)CLOCKS_PER_SEC);
// #endif
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