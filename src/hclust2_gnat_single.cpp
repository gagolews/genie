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


#include "hclust2_gnat_single.h"

using namespace Rcpp;
using namespace std;
using namespace boost;
using namespace DataStructures;


// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustGnatSingle::HClustGnatSingle(Distance* dist, RObject control) :
   opts(control), _root(NULL), _n(dist->getObjectCount()), _distance(dist),
   _indices(dist->getObjectCount()),
   neighborsCount(vector<size_t>(dist->getObjectCount(), 0)),
   minRadiuses(vector<double>(dist->getObjectCount(), -INFINITY)),
   // maxRadiuses(vector<double>(dist->getObjectCount(), INFINITY)),
   shouldFind(vector<bool>(dist->getObjectCount(), true)),
   nearestNeighbors(vector< deque<HeapNeighborItem> >(dist->getObjectCount())),
#ifdef GENERATE_STATS
   stats(HClustBiVpTreeStats()),
#endif
   ds(dist->getObjectCount())
{
#if VERBOSE > 5
   Rprintf("[%010.3f] building vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
#endif

   // starting _indices: random permutation of {0,1,...,_n-1}
   for(size_t i=0;i<_n;i++) _indices[i] = i;
   for(size_t i=_n-1; i>= 1; i--)
      swap(_indices[i], _indices[(size_t)(unif_rand()*(i+1))]);

   _root = buildFromPoints(opts.degree, 0, _n);
}


HClustGnatSingle::~HClustGnatSingle() {
// #if VERBOSE > 5
//       Rprintf("[%010.3f] destroying vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
// #endif
   if(_root) delete _root;
}


vector<size_t> HClustGnatSingle::chooseNewSplitPoints(size_t degree, size_t left, size_t right)
{
   Rcout <<"left= "<<left << " right= " << right << endl;
   const size_t candidatesTimes = opts.candidatesTimes;

   vector<size_t> splitPoints(degree);
   if(right-left <= degree)
   {
      for(size_t i=0;i<right-left;++i)
      {
         splitPoints[i] = _indices[left+i];
      }
      return splitPoints;
   }
   //wybierz losowo podzbior candidatesTimes * degree punktow
   size_t candidatesNumber = min(candidatesTimes*degree, right-left);
   Rcout <<"candidatesNumber= "<<candidatesNumber << endl;
   if(candidatesNumber > right-left)
   {
      for(size_t i=0;i<candidatesNumber;++i)
      {
         size_t tmpIndex = i + left + (size_t)(unif_rand()*(right-left-i));
         size_t tmp = _indices[tmpIndex];
         _indices[tmpIndex] = _indices[left+i];
         _indices[left+i] = tmp;
      }
   }
   //wybierz z tego podzbioru punkt
   size_t tmpIndex = left+(size_t)(unif_rand()*(candidatesNumber));
   size_t tmp = _indices[tmpIndex];
   _indices[tmpIndex] = _indices[left+0];
   _indices[left+0] = tmp;
   size_t lastAddedSplitPoint = splitPoints[0] = tmp;
   vector<vector<double>> dynamicProgrammingTable(degree-1); // we search for farest point only degree-1 times
   for(size_t i = 0; i<degree-1; ++i)
      dynamicProgrammingTable[i] = vector<double>(candidatesNumber);
   //znajdz jego najdalszego sasiada
   for(size_t i = 0; i<degree-1; ++i)
   {
      //znajdz kolejnego najdalszego do tych dwoch (minimum z dwoch dystansow musi byc najwieksze)
      //zrob  to z uzyciem programowania dynamicznego w O(nm), n-liczba kandydatow(candidatesTimes * degree), m-degree
      for(size_t j = 0; j<candidatesNumber; ++j)
      {
         if(j == lastAddedSplitPoint)
         {
            dynamicProgrammingTable[i][j] = 0;
            continue;
         }
         if(i > 0)
            dynamicProgrammingTable[i][j] = min(dynamicProgrammingTable[i-1][j], (*_distance)(lastAddedSplitPoint, _indices[left+j]));
         else
            dynamicProgrammingTable[i][j] = (*_distance)(lastAddedSplitPoint, _indices[left+j]);
      }
      size_t maxIndex = 0;
      size_t maxDist = dynamicProgrammingTable[i][0];
      for(size_t j = 1; j<candidatesNumber; ++j)
      {
         if(dynamicProgrammingTable[i][j] > maxDist)
         {
            maxIndex = j;
            maxDist = dynamicProgrammingTable[i][j];
         }
      }
      size_t tmp = _indices[left+maxIndex];
      _indices[left+maxIndex] = _indices[left+i];
      _indices[left+i] = tmp;
      lastAddedSplitPoint = splitPoints[i+1] = tmp;
   }

   return splitPoints;
}

HClustGnatSingleNode* HClustGnatSingle::createNonLeafNode(size_t degree, size_t left, size_t right,const vector<size_t>& splitPoints, const vector<size_t>& boundaries)
{
   HClustGnatSingleNode* node = new HClustGnatSingleNode();
   node->degree = degree;
   node->children = vector<HClustGnatSingleNode*>(degree);
   size_t childLeft = left+degree;
   size_t childRight;
   //assert: splitPoints.size() == degree
   //assert: boundaries.size() == degree //zawiera same prawe granice
   for(size_t i=0; i<degree; ++i)
   {
      childRight = boundaries[i];
      HClustGnatSingleNode* child = buildFromPoints(degree, childLeft, childRight);
      child->splitPointIndex = splitPoints[i];
      node->children[i] = child;
      childLeft = childRight;
   }

   return node;
}

vector<size_t> HClustGnatSingle::groupPointsToSplitPoints(const vector<size_t>& splitPoints, size_t left, size_t right)
{
   vector<vector<size_t>> groups(splitPoints.size());
   vector<size_t> boundaries(splitPoints.size());
   vector<double> distances(splitPoints.size());
   for(size_t i=left+splitPoints.size();i<right;++i)
   {
      size_t mySplitPointIndex = 0;
      double mySplitPointDistance = (*_distance)(splitPoints[0], _indices[i]);
      distances[0] = mySplitPointDistance;
      for(size_t j=1;j<splitPoints.size();++j)
      {
         double d = (*_distance)(splitPoints[j], _indices[i]);
         distances[j] = d;
         if(d < mySplitPointDistance)
         {
            mySplitPointIndex = j;
            mySplitPointDistance = d;
         }
      }
      groups[mySplitPointIndex].push_back(_indices[i]);

      for(size_t j=0;j<splitPoints.size();++j)
      {
         if(mySplitPointIndex == j)
            continue;
         auto rangeIterator = splitPointsRanges.find(SortedPoint(mySplitPointIndex, j));
         if(rangeIterator != splitPointsRanges.end())
         {
            if(distances[j] < rangeIterator->second.min)
            {
               rangeIterator->second.min = distances[j];
            }
            if(distances[j] > rangeIterator->second.max)
            {
               rangeIterator->second.max = distances[j];
            }
         }
         else
         {
            splitPointsRanges.emplace(SortedPoint(mySplitPointIndex, j), HClustGnatRange(distances[j], distances[j]));
         }
      }
   }
   size_t cumsum = 0;
   for(size_t i = 0; i<splitPoints.size(); ++i)
   {
      for(size_t j = 0; j<groups[i].size(); ++j)
      {
         _indices[left+splitPoints.size()+cumsum+j] = groups[i][j];
      }
      cumsum += groups[i].size();
      boundaries[i] = left+splitPoints.size()+cumsum;
   }
}

HClustGnatSingleNode* HClustGnatSingle::buildFromPoints(size_t degree, size_t left, size_t right)
{
#ifdef GENERATE_STATS
      ++stats.nodeCount;
#endif
   if(right - left <= opts.candidatesTimes*opts.degree) //@TODO: pomyslec, jaki tak naprawde mamy warunek stopu
   {
      Rcout << "tworze leaf" << endl;
      return new HClustGnatSingleNode(left, right);
   }
   Rcout << "split points" << endl;
   vector<size_t> splitPoints = chooseNewSplitPoints(degree, left, right);
   Rcout << "boundaries" << endl;
   vector<size_t> boundaries = groupPointsToSplitPoints(splitPoints, left, right); //@TODO: dobrze sie zastanowic, gdzie umieszczamy split pointy, aby nie szly w dol, gdzie sa granice!
   //@TODO: wybierac degree dziecka, zeby sie roznilo od degree aktualnego, jest w artykule
   Rcout << "tworze nonleaf" << endl;
   return createNonLeafNode(degree, left, right, splitPoints, boundaries);
}


void HClustGnatSingle::getNearestNeighborsFromMinRadiusRecursive( HClustGnatSingleNode* node, size_t index,
   size_t clusterIndex, double minR, double& maxR,
   std::priority_queue<HeapNeighborItem>& heap )
{
   // search within (minR, maxR]
   // if (node == NULL) return; // this should not happen
#ifdef GENERATE_STATS
   ++stats.nodeVisit;
#endif
   /*
   if (node->sameCluster) {
      if (node->splitPointIndex == SIZE_MAX) {
         if (ds.find_set(_indices[node->left]) == clusterIndex) return;
      } else {
         if (ds.find_set(node->splitPointIndex) == clusterIndex) return;
      }
   }*/

   if (node->degree == SIZE_MAX) // leaf
   {
      /*if (node->sameCluster)
      {
         for (size_t i=node->left; i<node->right; i++)
         {
            if (index >= _indices[i]) continue;
            double dist2 = (*_distance)(index, _indices[i]);
            if (dist2 > maxR || dist2 <= minR) continue;

            if (heap.size() >= opts.maxNNPrefetch) {
               if (dist2 < maxR) {
                  while (!heap.empty() && heap.top().dist == maxR) {
                     heap.pop();
                  }
               }
            }
            heap.push( HeapNeighborItem(_indices[i], dist2) );
            maxR = heap.top().dist;
         }
      }
      else*/
      {
         //size_t commonCluster = ds.find_set(_indices[node->left]);
         for (size_t i=node->left; i<node->right; i++)
         {
            //size_t currentCluster = ds.find_set(_indices[i]);
            //if (currentCluster != commonCluster) commonCluster = SIZE_MAX;
            //if (currentCluster == clusterIndex) continue;

            if (index >= _indices[i]) continue;

            double dist2 = (*_distance)(index, _indices[i]);
            if (dist2 > maxR || dist2 <= minR) continue;

            if (heap.size() >= opts.maxNNPrefetch) {
               if (dist2 < maxR) {
                  while (!heap.empty() && heap.top().dist == maxR) {
                     heap.pop();
                  }
               }
            }
            heap.push( HeapNeighborItem(_indices[i], dist2) );
            maxR = heap.top().dist;
         }
         //if (commonCluster != SIZE_MAX) node->sameCluster = true;
      }
      return;
   }
   // else // not a leaf
   //1. z artykulu
   vector<bool> shouldVisit(node->degree, true);
   for(size_t i=0;i<node->degree; ++i) //4. z artykulu
   {
      if(shouldVisit[i])
      {
         size_t pi = node->children[i]->splitPointIndex;
         //2. z artykulu
         double dist = (*_distance)(node->children[i]->splitPointIndex, index);

         if (ds.find_set(node->splitPointIndex) != clusterIndex && index < node->splitPointIndex) {
               if (dist <= maxR && dist > minR) {
                  if (heap.size() >= opts.maxNNPrefetch && dist < maxR) {
                     while (!heap.empty() && heap.top().dist == maxR) {
                        heap.pop();
                     }
                  }
                  heap.push( HeapNeighborItem(node->splitPointIndex, dist) );
                  maxR = heap.top().dist;
               }
            }
         //3. z artykulu
         for(size_t j=0;j<node->degree;++j)
         {
            if(i != j && shouldVisit[j])
            {
               size_t pj = node->children[j]->splitPointIndex;
               auto rangeIterator = splitPointsRanges.find(SortedPoint(pi, pj));
               //assert: rangeIterator != splitPointsRanges.end()
               HClustGnatRange range = rangeIterator->second;
               double leftRange = dist-maxR;
               double rightRange = dist+maxR;
               if(leftRange <= range.max && range.min <= rightRange) //http://world.std.com/~swmcd/steven/tech/interval.html
               {//they intersect
                  ;
               }
               else
               {//disjoint
                  shouldVisit[j] = false;
               }
            }
         }

      }
   }
   //5. z artykulu
   for(size_t i=0;i<node->degree; ++i) //4. z artykulu
   {
      if(shouldVisit[i])
      {
         getNearestNeighborsFromMinRadiusRecursive(node->children[i], index, clusterIndex, minR, maxR, heap);
      }
   }

   //@TODO: robic same cluster
   /*if (   !node->sameCluster
      && (!node->ll || node->ll->sameCluster)
      && (!node->rl || node->rl->sameCluster)
#ifndef USE_ONEWAY_VPTREE
      && (!node->lr || node->lr->sameCluster)
      && (!node->rr || node->rr->sameCluster)
#endif
      )
   {
      size_t commonCluster = SIZE_MAX;
      if (node->ll) {
         size_t currentCluster = ds.find_set((node->ll->splitPointIndex == SIZE_MAX)?_indices[node->ll->left]:node->ll->splitPointIndex);
         if (commonCluster == SIZE_MAX) commonCluster = currentCluster;
         else if (currentCluster != commonCluster) return;
      }
      if (node->rl) {
         size_t currentCluster = ds.find_set((node->rl->splitPointIndex == SIZE_MAX)?_indices[node->rl->left]:node->rl->splitPointIndex);
         if (commonCluster == SIZE_MAX) commonCluster = currentCluster;
         else if (currentCluster != commonCluster) return;
      }
#ifndef USE_ONEWAY_VPTREE
      if (node->lr) {
         size_t currentCluster = ds.find_set((node->lr->splitPointIndex == SIZE_MAX)?_indices[node->lr->left]:node->lr->splitPointIndex);
         if (commonCluster == SIZE_MAX) commonCluster = currentCluster;
         else if (currentCluster != commonCluster) return;
      }
      if (node->rr) {
         size_t currentCluster = ds.find_set((node->rr->splitPointIndex == SIZE_MAX)?_indices[node->rr->left]:node->rr->splitPointIndex);
         if (commonCluster == SIZE_MAX) commonCluster = currentCluster;
         else if (currentCluster != commonCluster) return;
      }
#endif
      node->sameCluster = true;
   }
   */
}

void HClustGnatSingle::FindNeighborTest(size_t index, double R)
{
   std::priority_queue<HeapNeighborItem> heap;
   size_t clusterIndex = ds.find_set(index);
   double _tau = R;
   getNearestNeighborsFromMinRadiusRecursive( _root, index, clusterIndex, minRadiuses[index], _tau, heap );
   while( !heap.empty() ) {
      Rcout << heap.top().index << ", " << heap.top().dist << endl;
      heap.pop();
   }
}


HeapNeighborItem HClustGnatSingle::getNearestNeighbor(size_t index)
{
#if VERBOSE > 5
   // Rprintf(".");
#endif
//   This is not faster:
//       while(!nearestNeighbors[index].empty())
//       {
//          size_t sx = ds.find_set(index);
//          auto res = nearestNeighbors[index].front();
//          size_t sy = ds.find_set(res.index);
//          nearestNeighbors[index].pop_front();
//          if (sx != sy) {
//             return res;
//          }
//          // else just go on removing items
//       }


   if(shouldFind[index] && nearestNeighbors[index].empty())
   {
      std::priority_queue<HeapNeighborItem> heap;
      size_t clusterIndex = ds.find_set(index);

      double _tau = INFINITY;//maxRadiuses[index];

//       THIS IS SLOWER:
//          double _tau = (*_distance)(index,
//             ds.getClusterNext(clusterIndex))
//          );

//       THIS IS SLOWER TOO:
//          size_t test = (size_t)(index+unif_rand()*(_n-index));
//          if (ds.find_set(test) != clusterIndex)
//             _tau = (*_distance)(index, test);

#ifdef GENERATE_STATS
      ++stats.nnCals;
#endif
      getNearestNeighborsFromMinRadiusRecursive( _root, index, clusterIndex, minRadiuses[index], _tau, heap );
      while( !heap.empty() ) {
         nearestNeighbors[index].push_front(heap.top());
         heap.pop();
      }
      // maxRadiuses[index] = INFINITY;
      size_t newNeighborsCount = nearestNeighbors[index].size();

      neighborsCount[index] += newNeighborsCount;
      if(neighborsCount[index] > _n - index || newNeighborsCount == 0)
         shouldFind[index] = false;

      if(newNeighborsCount > 0)
         minRadiuses[index] = nearestNeighbors[index].back().dist;
   }

   if(!nearestNeighbors[index].empty())
   {
#ifdef GENERATE_STATS
      ++stats.nnCount;
#endif
      auto res = nearestNeighbors[index].front();
      nearestNeighbors[index].pop_front();
      return res;
   }
   else
   {
      return HeapNeighborItem(SIZE_MAX,-INFINITY);
      //stop("nie ma sasiadow!");
   }
}


NumericMatrix HClustGnatSingle::compute()
{
   NumericMatrix ret(_n-1, 2);
   priority_queue<HeapHierarchicalItem> pq;

   // INIT: Pre-fetch a few nearest neighbors for each point
#if VERBOSE > 5
   Rprintf("[%010.3f] prefetching NNs\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   for (size_t i=0; i<_n; i++)
   {
#if VERBOSE > 7
      if (i % 1024 == 0) Rprintf("\r             prefetch NN: %d/%d", i, _n-1);
#endif
      Rcpp::checkUserInterrupt(); // may throw an exception, fast op
      HeapNeighborItem hi=getNearestNeighbor(i);

      if (hi.index != SIZE_MAX)
      {
         pq.push(HeapHierarchicalItem(i, hi.index, hi.dist));
      }
   }
#if VERBOSE > 7
   Rprintf("\r             prefetch NN: %d/%d\n", _n-1, _n-1);
#endif
#if VERBOSE > 5
   Rprintf("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);
#endif

   size_t i = 0;
   while(i < _n - 1)
   {
      //Rcout << "iteracja " << i << endl;
      //Rcout << "pq size = " << pq.size()<< endl;
      HeapHierarchicalItem hhi = pq.top();
      pq.pop();

      size_t s1 = ds.find_set(hhi.index1);
      size_t s2 = ds.find_set(hhi.index2);
      if (s1 != s2)
      {
         Rcpp::checkUserInterrupt(); // may throw an exception, fast op

         ret(i,0)=(double)hhi.index1;
         ret(i,1)=(double)hhi.index2;
         ++i;
         ds.link(s1, s2);
      }
#if VERBOSE > 7
      if (i % 1024 == 0) Rprintf("\r             %d / %d", i+1, _n);
#endif

      // ASSERT: hhi.index1 < hhi.index2
      HeapNeighborItem hi=getNearestNeighbor(hhi.index1);
      if(hi.index != SIZE_MAX)
         pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));
   }
#if VERBOSE > 7
   Rprintf("\r             %d / %d\n", _n, _n);
#endif
#if VERBOSE > 5
   Rprintf("[%010.3f] generating output matrix\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   Rcpp::checkUserInterrupt();

   MergeMatrixGenerator mmg(ret.nrow());
   return mmg.generateMergeMatrix(ret);
}




// [[Rcpp::export]]
void hclust_gnat_single_test(RObject distance, RObject objects, RObject control=R_NilValue, int index=0, double R=0) {
   RObject result(R_NilValue);
   DataStructures::Distance* dist = DataStructures::Distance::createDistance(distance, objects);
   DataStructures::HClustGnatSingle hclust(dist, control);
   //RObject merge = hclust.compute();
   Rcpp::stop("GNAT zbudowany!");
   hclust.FindNeighborTest(index, R);
   if (dist) delete dist;
}

// [[Rcpp::export]]
RObject hclust_gnat_single(RObject distance, RObject objects, RObject control=R_NilValue) {
#if VERBOSE > 5
   Rprintf("[%010.3f] starting timer\n", clock()/(double)CLOCKS_PER_SEC);
#endif
   RObject result(R_NilValue);
   DataStructures::Distance* dist = DataStructures::Distance::createDistance(distance, objects);

   try {
      /* Rcpp::checkUserInterrupt(); may throw an exception */
      DataStructures::HClustGnatSingle hclust(dist, control);
      RObject merge = hclust.compute();
      result = Rcpp::as<RObject>(List::create(
         _["merge"]  = merge,
         _["height"] = R_NilValue,
         _["order"]  = R_NilValue,
         _["labels"] = R_NilValue,
         _["call"]   = R_NilValue,
         _["method"] = "single",
         _["dist.method"] = R_NilValue,
         _["stats"] = List::create(
            _["vptree"] = hclust.getStats().toR(),
            _["distance"] = dist->getStats().toR()
         ),
         _["control"] = List::create(
            _["vptree"] = hclust.getOptions().toR()
         )
      ));
      result.attr("class") = "hclust";
      //hclust.print();
   }
   catch(...) {

   }

   if (dist) delete dist;
#if VERBOSE > 5
   Rprintf("[%010.3f] done\n", clock()/(double)CLOCKS_PER_SEC);
#endif
   if (Rf_isNull(result)) stop("stopping on error or explicit user interrupt");
   return result;
}


