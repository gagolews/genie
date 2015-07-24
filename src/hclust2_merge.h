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



#include <vector>
#include <cstdint>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

namespace DataStructures{
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
       MergeMatrixGenerator(size_t n);
       NumericMatrix generateMergeMatrix(const NumericMatrix& x);
    private:
       void initialize(const NumericMatrix& x, size_t k, size_t clusterNumber);
       void writeRowIfSingleElements(NumericMatrix& y, size_t k);
       void writeRowIfClusterElements(NumericMatrix& y, size_t k);
       size_t findMyParent(size_t s, size_t clusterNumber);
    };


 }
