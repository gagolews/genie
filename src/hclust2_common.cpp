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

#include "hclust2_common.h"

using namespace grup;

HClustOptions* NNHeap::opts = NULL;


HClustOptions::HClustOptions(Rcpp::RObject control) {
   maxLeavesElems = DEFAULT_MAX_LEAVES_ELEMS;
   maxNNPrefetch = DEFAULT_MAX_NN_PREFETCH;
   maxNNMerge = DEFAULT_MAX_NN_MERGE;
   minNNPrefetch = DEFAULT_MIN_NN_PREFETCH;
   minNNMerge = DEFAULT_MIN_NN_MERGE;
   vpSelectScheme = DEFAULT_VP_SELECT_SCHEME;
   vpSelectCand = DEFAULT_VP_SELECT_CAND;
   vpSelectTest = DEFAULT_VP_SELECT_TEST;
   nodesVisitedLimit = DEFAULT_NODES_VISITED_LIMIT;
   thresholdGini = DEFAULT_THRESHOLD_GINI;
   useVpTree = DEFAULT_USEVPTREE;
   useMST = DEFAULT_USEMST;

   if (!Rf_isNull((SEXP)control)) {
      Rcpp::List control2(control);

      // if (control2.containsElementNamed("exemplar")) {
      //    Rcpp::CharacterVector ex = Rcpp::clone(Rcpp::as<Rcpp::CharacterVector>(control2["exemplar"]));
      //    const char* exc = CHAR(STRING_ELT((SEXP)ex, 0));
      //    exemplar = std::string(exc);
      //    //Rcpp::CharacterVector exemplar = Rcpp::CharacterVector(control2["exemplar"]);
      // }

      if (control2.containsElementNamed("maxLeavesElems")) {
         maxLeavesElems = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["maxLeavesElems"])[0];
      }

      if (control2.containsElementNamed("maxNNPrefetch")) {
         maxNNPrefetch = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["maxNNPrefetch"])[0];
      }

      if (control2.containsElementNamed("maxNNMerge")) {
         maxNNMerge = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["maxNNMerge"])[0];
      }

      if (control2.containsElementNamed("minNNPrefetch")) {
         minNNPrefetch = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["minNNPrefetch"])[0];
      }

      if (control2.containsElementNamed("minNNMerge")) {
         minNNMerge = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["minNNMerge"])[0];
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

      if (control2.containsElementNamed("nodesVisitedLimit")) {
         nodesVisitedLimit = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["nodesVisitedLimit"])[0];
      }

      if (control2.containsElementNamed("thresholdGini")) {
         thresholdGini = (double)Rcpp::as<Rcpp::NumericVector>(control2["thresholdGini"])[0];
      }

      if (control2.containsElementNamed("useVpTree")) {
         useVpTree = (bool)Rcpp::as<Rcpp::LogicalVector>(control2["useVpTree"])[0];
      }

      if (control2.containsElementNamed("useMST")) {
         useMST = (bool)Rcpp::as<Rcpp::LogicalVector>(control2["useMST"])[0];
      }
   }

   if (thresholdGini < 0.0 || thresholdGini > 1.0) {
      thresholdGini = DEFAULT_THRESHOLD_GINI;
      Rf_warning("wrong thresholdGini value. using default");
   }
   if (maxLeavesElems < 2 || maxLeavesElems > 512) {
      maxLeavesElems = DEFAULT_MAX_LEAVES_ELEMS;
      Rf_warning("wrong maxLeavesElems value. using default");
   }
   if (maxNNPrefetch < 1) {
      maxNNPrefetch = DEFAULT_MAX_NN_PREFETCH;
      Rf_warning("wrong maxNNPrefetch value. using default");
   }
   if (maxNNMerge < 1) {
      STOPIFNOT(maxNNPrefetch >= 1 && maxNNPrefetch != SIZE_MAX);
      maxNNMerge = DEFAULT_MAX_NN_MERGE;
      Rf_warning("wrong maxNNMerge value. using default");
   }
   if (minNNPrefetch < 1) {
      minNNPrefetch = DEFAULT_MIN_NN_PREFETCH;
      Rf_warning("wrong minNNPrefetch value. using default");
   }
   if (minNNMerge < 1) {
      STOPIFNOT(minNNPrefetch >= 1 && minNNPrefetch != SIZE_MAX);
      minNNMerge = DEFAULT_MIN_NN_MERGE;
      Rf_warning("wrong minNNMerge value. using default");
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
   if (nodesVisitedLimit < 1) {
      nodesVisitedLimit = DEFAULT_NODES_VISITED_LIMIT;
      Rf_warning("wrong nodesVisitedLimit value. using default");
   }
}


Rcpp::NumericVector HClustOptions::toR() const
{
   return Rcpp::NumericVector::create(
      Rcpp::_["maxLeavesElems"]     = maxLeavesElems,
      Rcpp::_["maxNNPrefetch"]      = maxNNPrefetch,
      Rcpp::_["maxNNMerge"]         = maxNNMerge,
      Rcpp::_["minNNPrefetch"]      = minNNPrefetch,
      Rcpp::_["minNNMerge"]         = minNNMerge,
      Rcpp::_["vpSelectScheme"]     = vpSelectScheme,
      Rcpp::_["vpSelectCand"]       = vpSelectCand,
      Rcpp::_["vpSelectTest"]       = vpSelectTest,
      Rcpp::_["nodesVisitedLimit"]  = nodesVisitedLimit,
      Rcpp::_["thresholdGini"]      = thresholdGini,
      Rcpp::_["useVpTree"]          = useVpTree,
      Rcpp::_["useMST"]             = useMST
   );
}


HClustStats::HClustStats() :
   nodeCount(0), leafCount(0), nodeVisit(0), nnCals(0), nnCount(0),
   medoidOldNew(0), medoidUpdateCount(0) {}

HClustStats::~HClustStats() {
   #if VERBOSE > 0
   Rprintf("             vp-tree: nodeCount=%.0f, leafCount=%.0f, nodeVisit=%.0f, nnCals=%.0f, nnCount=%.0f, medoidUpdateCount=%.0f, medoidOldNew=%.0f\n",
      (double)nodeCount, (double)leafCount, (double)nodeVisit, (double)nnCals, (double)nnCount,(double)medoidUpdateCount, (double)medoidOldNew);
   #endif
}


Rcpp::NumericVector HClustStats::toR() const {
   return Rcpp::NumericVector::create(
      Rcpp::_["nodeCount"]
         = (nodeCount>0)?(double)nodeCount:NA_REAL,
      Rcpp::_["leafCount"]
         = (leafCount>0)?(double)leafCount:NA_REAL,
      Rcpp::_["nodeVisit"]
         = (nodeVisit>0)?(double)nodeVisit:NA_REAL,
      Rcpp::_["nnCals"]
         = (nnCals>0)?(double)nnCals:NA_REAL,
      Rcpp::_["nnCount"]
         = (nnCount>0)?(double)nnCount:NA_REAL,
      Rcpp::_["medoidUpdateCount"]
         = (medoidUpdateCount>0)?medoidUpdateCount:NA_REAL,
      Rcpp::_["medoidOldNew"]
         = (medoidOldNew>0)?medoidOldNew:NA_REAL
   );
}
