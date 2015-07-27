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

#include "hclust2_common.h"

using namespace DataStructures;


HClustBiVpTreeOptions::HClustBiVpTreeOptions(Rcpp::RObject control) {
   maxLeavesElems = DEFAULT_MAX_LEAVES_ELEMS;
   maxNNPrefetch  = DEFAULT_MAX_NN_PREFETCH;
   vpSelectScheme = DEFAULT_VP_SELECT_SCHEME;
   vpSelectCand   = DEFAULT_VP_SELECT_CAND;
   vpSelectTest   = DEFAULT_VP_SELECT_TEST;

   if (!Rf_isNull((SEXP)control)) {
      Rcpp::List control2(control);

      if (control2.containsElementNamed("maxLeavesElems")) {
         maxLeavesElems = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["maxLeavesElems"])[0];
      }

      if (control2.containsElementNamed("maxNNPrefetch")) {
         maxNNPrefetch = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["maxNNPrefetch"])[0];
      }

      if (control2.containsElementNamed("vpSelectScheme")) {
         vpSelectScheme = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["vpSelectScheme"])[0];
      }

      if (control2.containsElementNamed("vpSelectCand")) {
         vpSelectCand = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["vpSelectCand"])[0];
      }

      if (control2.containsElementNamed("vpSelectTest")) {
         vpSelectTest = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["vpSelectTest"])[0];
      }
   }

   if (maxLeavesElems < 2 || maxLeavesElems > 64) {
      maxLeavesElems = DEFAULT_MAX_LEAVES_ELEMS;
      Rf_warning("wrong maxLeavesElems value. using default");
   }
   if (maxNNPrefetch < 1 || maxNNPrefetch > 16) {
      maxNNPrefetch = DEFAULT_MAX_NN_PREFETCH;
      Rf_warning("wrong maxNNPrefetch value. using default");
   }
   if (vpSelectScheme < 1 || vpSelectScheme > 3) {
      vpSelectScheme = DEFAULT_VP_SELECT_SCHEME;
      Rf_warning("wrong vpSelectScheme value. using default");
   }
   if (vpSelectCand < 1 || vpSelectCand > 128) {
      vpSelectCand = DEFAULT_VP_SELECT_CAND;
      Rf_warning("wrong vpSelectCand value. using default");
   }
   if (vpSelectTest < 1 || vpSelectTest > 128) {
      vpSelectTest = DEFAULT_VP_SELECT_TEST;
      Rf_warning("wrong vpSelectTest value. using default");
   }
}


Rcpp::NumericVector HClustBiVpTreeOptions::toR() const
{
   return Rcpp::NumericVector::create(
      Rcpp::_["maxLeavesElems"] = maxLeavesElems,
      Rcpp::_["maxNNPrefetch"] = maxNNPrefetch,
      Rcpp::_["vpSelectScheme"] = vpSelectScheme,
      Rcpp::_["vpSelectCand"] = vpSelectCand,
      Rcpp::_["vpSelectTest"] = vpSelectTest
   );
}


HClustBiVpTreeStats::HClustBiVpTreeStats() :
   nodeCount(0), nodeVisit(0), nnCals(0), nnCount(0) {}

HClustBiVpTreeStats::~HClustBiVpTreeStats() {
   #if VERBOSE > 0
   Rprintf("             vp-tree: nodeCount=%.0f, nodeVisit=%.0f, nnCals=%.0f, nnCount=%.0f\n",
      (double)nodeCount, (double)nodeVisit, (double)nnCals, (double)nnCount);
   #endif
}


Rcpp::NumericVector HClustBiVpTreeStats::toR() const {
   return Rcpp::NumericVector::create(
      Rcpp::_["nodeCount"]
         = (nodeCount>0)?(double)nodeCount:NA_REAL,
      Rcpp::_["nodeVisit"]
         = (nodeVisit>0)?(double)nodeVisit:NA_REAL,
      Rcpp::_["nnCals"]
         = (nnCals>0)?(double)nnCals:NA_REAL,
      Rcpp::_["nnCount"]
         = (nnCount>0)?(double)nnCount:NA_REAL
   );
}
