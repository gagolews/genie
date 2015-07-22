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
