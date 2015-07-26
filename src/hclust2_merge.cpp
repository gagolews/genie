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


#include "hclust2_merge.h"

using namespace DataStructures;

MergeMatrixGenerator::MergeMatrixGenerator(size_t n)
    : n(n), elements(n+2,0), parents(n+2,0),
      i(SIZE_MAX), j(SIZE_MAX), si(SIZE_MAX), sj(SIZE_MAX)
{
}

NumericMatrix MergeMatrixGenerator::generateMergeMatrix(const NumericMatrix& x)
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

void MergeMatrixGenerator::initialize(const NumericMatrix& x, size_t k, size_t clusterNumber)
{
 if(x(k,0) < 0 || x(k,1) < 0) stop("There are negative values. You probably gave merged matrix.");
 i = (size_t)x(k,0)+1;
 j = (size_t)x(k,1)+1;
 si = elements[i];
 sj = elements[j];
 elements[i] = clusterNumber;
 elements[j] = clusterNumber;
}
void MergeMatrixGenerator::writeRowIfSingleElements(NumericMatrix& y, size_t k)
{
 if(si == 0)
    y(k,0)=-(double)i;

 if(sj == 0)
    y(k,1)=-(double)j;
}

void MergeMatrixGenerator::writeRowIfClusterElements(NumericMatrix& y, size_t k)
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

size_t MergeMatrixGenerator::findMyParent(size_t s, size_t clusterNumber)
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


// [[Rcpp::export]]
NumericMatrix generateMergeMatrix(NumericMatrix x) {
   DataStructures::MergeMatrixGenerator mmg(x.nrow());
   return mmg.generateMergeMatrix(x);
}