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

#ifndef __DEFS_H
#define __DEFS_H

#include <cstdint> /* SIZE_MAX, C++11 */
#include <limits>  /* INFINITY, C++11 */
#include <cmath>   /* INFINITY, C++11 */

#ifndef INFINITY
#define INFINITY (std::numeric_limits<double>::infinity())
#endif

#ifndef SIZE_MAX
#define SIZE_MAX (std::numeric_limits<size_t>::max())
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _OPENMP
#define MASTER_OR_SINGLE_THREAD (omp_get_thread_num() == 0)
#else
#define MASTER_OR_SINGLE_THREAD (1)
#endif

#ifdef _OPENMP
#define OPENMP_ONLY(x) {x;}
#else
#define OPENMP_ONLY(x)  ;
#endif


// ---------------------------------------------------------------------
// #define HASHMAP_ENABLED



#if 0
/* DEBUG MODE */
#define VERBOSE 1
#define GENERATE_STATS
#define DISJOINT_SETS_DEBUG
#define MEASURE_MEM_USE
#else
/* PRODUCTION USE */
#define VERBOSE 0
#endif


#define DEFAULT_MAX_LEAVES_ELEMS 4
#define DEFAULT_MAX_NN_PREFETCH 256
#define DEFAULT_MAX_NN_MERGE maxNNPrefetch
#define DEFAULT_MIN_NN_PREFETCH 20
#define DEFAULT_MIN_NN_MERGE minNNPrefetch
#define DEFAULT_VP_SELECT_SCHEME 3
#define DEFAULT_VP_SELECT_CAND 5
#define DEFAULT_VP_SELECT_TEST 12
#define DEFAULT_NODES_VISITED_LIMIT SIZE_MAX
#define DEFAULT_THRESHOLD_GINI 0.3
#define DEFAULT_USEVPTREE false
#define DEFAULT_USEMST true
// #define DEFAULT_GNAT_DEGREE 50
// #define DEFAULT_GNAT_CANDIDATES_TIMES 3
// #define DEFAULT_GNAT_MIN_DEGREE 2
// #define DEFAULT_GNAT_MAX_DEGREE 200
// #define DEFAULT_GNAT_MAX_TIMES_DEGREE 5
// #define DEFAULT_EXEMPLAR_UPDATE_METHOD 2
// #define DEFAULT_EXEMPLAR_MAX_LEAVES_ELEMS 32
// #define DEFAULT_IS_CURSE_OF_DIMENSIONALITY false

// ---------------------------------------------------------------------------

#if VERBOSE > 0
#define STOPIFNOT(EXPR) { if (!(EXPR)) Rprintf("\a*** Assert failed: " #EXPR " at %s, line %d ***\n", __FILE__, __LINE__); }
#else
#define STOPIFNOT(EXPR) { }
#endif

// ---------------------------------------------------------------------------


#if VERBOSE > 0
#define RCOUT(msg, verlvl) if ((verlvl) <= VERBOSE) Rcout << msg << endl;
#else
#define RCOUT(msg, verlvl) { }
#endif

/*
 * example use:
 *  int a = 5;
 *  RCOUT("a=" << a, 0)
 *
 */
// ---------------------------------------------------------------------


/* TO DO: can we do it more intelligently? */

#if VERBOSE < 1
#define MESSAGE_1(...)  { }
#else
#define MESSAGE_1(...) Rprintf(__VA_ARGS__)
#endif

#if VERBOSE < 2
#define MESSAGE_2(...)  { }
#else
#define MESSAGE_2(...) Rprintf(__VA_ARGS__)
#endif

#if VERBOSE < 3
#define MESSAGE_3(...) { }
#else
#define MESSAGE_3(...) Rprintf(__VA_ARGS__)
#endif

#if VERBOSE < 4
#define MESSAGE_4(...) { }
#else
#define MESSAGE_4(...) Rprintf(__VA_ARGS__)
#endif

#if VERBOSE < 5
#define MESSAGE_5(...) { }
#else
#define MESSAGE_5(...) Rprintf(__VA_ARGS__)
#endif

#if VERBOSE < 6
#define MESSAGE_6(...) { }
#else
#define MESSAGE_6(...) Rprintf(__VA_ARGS__)
#endif

#if VERBOSE < 7
#define MESSAGE_7(...) { }
#else
#define MESSAGE_7(...) Rprintf(__VA_ARGS__)
#endif

#if VERBOSE < 8
#define MESSAGE_8(...) { }
#else
#define MESSAGE_8(...) Rprintf(__VA_ARGS__)
#endif

#if VERBOSE < 9
#define MESSAGE_9(...) { }
#else
#define MESSAGE_9(...) Rprintf(__VA_ARGS__)
#endif

#endif
