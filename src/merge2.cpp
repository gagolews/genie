/*
 * merge2.cpp
 *
 *  Created on: 21 lip 2015
 *      Author: bartoszukm
 */

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
#include <deque>
#include <exception>
#include <unordered_map>
#include <string>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/functional/hash.hpp>

using namespace Rcpp;
using namespace std;
using namespace boost;

namespace experiments
{

class MergeMatrixGenerator
   {
      size_t n;
      vector<size_t> elements;
      vector<size_t> parents;

      size_t i;
      size_t j;
      size_t si;
      size_t sj;
   public:
      MergeMatrixGenerator(size_t n)
         : n(n), elements(n+2,0), parents(n+2,0),
           i(SIZE_MAX), j(SIZE_MAX), si(SIZE_MAX), sj(SIZE_MAX)
         {
         }

      NumericMatrix generateMergeMatrix(const NumericMatrix& x)
      {
         // x has 0-based indices
         if (x.ncol() != 2) stop("x should have 2 columns");
         NumericMatrix y(n, 2);

         size_t clusterNumber = 1;
         for (size_t k=0; k<n; ++k, ++clusterNumber)
         {
            initialize(x, k, clusterNumber);
            writeRowIfSingleElements(y,k);
            si = findMyParent(si, clusterNumber);
            sj = findMyParent(sj, clusterNumber);
            writeRowIfClusterElements(y,k);
         }

         return y;
      }
   private:
      void initialize(const NumericMatrix& x, size_t k, size_t clusterNumber)
      {
         if(x(k,0) < 0 || x(k,1) < 0) stop("There are negative values. You probably gave merged matrix.");
         i = (size_t)x(k,0)+1;
         j = (size_t)x(k,1)+1;

         si = elements[i];
         sj = elements[j];
         elements[i] = clusterNumber;
         elements[j] = clusterNumber;
      }
      void writeRowIfSingleElements(NumericMatrix& y, size_t k)
      {
         if(si == 0)
            y(k,0)=-(double)i;

         if(sj == 0)
            y(k,1)=-(double)j;
      }

      void writeRowIfClusterElements(NumericMatrix& y, size_t k)
      {
         if(y(k,0) == 0)
         {
            y(k,0) = (double)si;
         }

         if(y(k,1) == 0)
         {
            y(k,1) = (double)sj;
         }
      }

      size_t findMyParent(size_t s, size_t clusterNumber)
      {
         while(parents[s] != 0)
         {
            size_t sinew = parents[s];
            parents[s] = clusterNumber;
            s = sinew;
         }
         if(s != 0) parents[s] = clusterNumber;
         return s;
      }
   };
}

// [[Rcpp::export]]
NumericMatrix generateMergeMatrix2(NumericMatrix x) {
   experiments::MergeMatrixGenerator mmg(x.nrow());
   return mmg.generateMergeMatrix(x);
}

