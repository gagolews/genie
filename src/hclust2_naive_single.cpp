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


#include "hclust2_naive_single.h"


using namespace Rcpp;
using namespace std;
using namespace boost;
using namespace DataStructures;


// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustNaiveSingle::HClustNaiveSingle(Distance* dist, RObject control) :
   HClustNNbasedSingle(dist, control)
{

}


HClustNaiveSingle::~HClustNaiveSingle()
{

}


void HClustNaiveSingle::getNearestNeighborsFromMinRadius(size_t index,
      size_t clusterIndex, double minR, NNHeap& nnheap)
{
   double maxR = INFINITY;
   for (size_t i=index+1; i<_n; ++i) {
      size_t currentCluster = ds.find_set(i);
      if (currentCluster == clusterIndex) continue;
      double dist2 = (*_distance)(_indices[index], _indices[i]); // the slow part
      if (dist2 > maxR || dist2 <= minR) continue;
      nnheap.insert(i, dist2, maxR);
   }
}
