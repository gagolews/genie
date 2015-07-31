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

#define VERBOSE 0

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
   printIndices();

   _root = buildFromPoints(opts.degree, 0, _n);
}


HClustGnatSingle::~HClustGnatSingle() {
// #if VERBOSE > 5
//       Rprintf("[%010.3f] destroying vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
// #endif
   if(_root) delete _root;
}

void HClustGnatSingle::printIndices()
{
#if VERBOSE > 11
   for(size_t i=0;i<_indices.size();++i)
      Rcout << _indices[i] << ", ";
   Rcout << endl;
#endif
   ;
}

vector<size_t> HClustGnatSingle::chooseNewSplitPoints(size_t degree, size_t left, size_t right)
{
   RCOUT("left= "<<left << " right= " << right,11)
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
   RCOUT("candidatesNumber= "<<candidatesNumber,11)
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
   size_t tmpIndex = (size_t)(unif_rand()*(candidatesNumber));
   //size_t tmp = _indices[tmpIndex];
   //_indices[tmpIndex] = _indices[left+0];
   //_indices[left+0] = tmp;
   RCOUT("tmpIndex = " << tmpIndex,11)
   size_t lastAddedSplitPoint = splitPoints[0] = tmpIndex;
   RCOUT("indices[tmpIndex] = " << _indices[left+tmpIndex],11)
   vector<vector<double>> dynamicProgrammingTable(degree-1); // we search for farest point only degree-1 times
   for(size_t i = 0; i<degree-1; ++i)
      dynamicProgrammingTable[i] = vector<double>(candidatesNumber);
   unordered_map<SortedPoint, HClustGnatRange> splitPointsRanges_small;
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
         RCOUT("licze odleglosc miedzy " << _indices[left+lastAddedSplitPoint]<< " a " << _indices[left+j], 11)
         double d = (*_distance)(_indices[left+lastAddedSplitPoint], _indices[left+j]);
         splitPointsRanges_small.emplace(SortedPoint(_indices[left+lastAddedSplitPoint], _indices[left+j]), HClustGnatRange(d,d));
         if(i > 0)
         {
            dynamicProgrammingTable[i][j] = min(dynamicProgrammingTable[i-1][j], d);
         }
         else
         {
            dynamicProgrammingTable[i][j] = d;
         }
      }
      size_t maxIndex = 0;
      double maxDist = dynamicProgrammingTable[i][0];
      //Rcout <<"dynamicProgrammingTable["<<i<<"][0]=" <<maxDist << endl;
      for(size_t j = 1; j<candidatesNumber; ++j)
      {
         //Rcout <<"dynamicProgrammingTable["<<i<<"]["<<j<<"]=" <<dynamicProgrammingTable[i][j] << endl;
         if(dynamicProgrammingTable[i][j] > maxDist)
         {
            maxIndex = j;
            maxDist = dynamicProgrammingTable[i][j];
            //Rcout << "chosen " <<_indices[left+maxIndex] << " i.e. " << maxIndex<< "with distance" <<maxDist << endl;
         }
      }
      //size_t tmp = _indices[left+maxIndex];
      //_indices[left+maxIndex] = _indices[left+i];
      //_indices[left+i] = tmp;
      lastAddedSplitPoint = splitPoints[i+1] = maxIndex;
   }

#if VERBOSE > 11
   Rcout << "gole indeksy" << endl;
   for(size_t i = 0; i < splitPoints.size(); ++i)
   {
      Rcout << splitPoints[i]  <<" ";
   }
   Rcout << endl;

   Rcout << "indeksy punktow" << endl;
   for(size_t i = 0; i < splitPoints.size(); ++i)
   {
      Rcout << _indices[left+splitPoints[i]]  <<" ";
   }
   Rcout << endl;
#endif

   vector<bool> canSwap(splitPoints.size(), true);
   sort(splitPoints.begin(), splitPoints.end());
   int index = 0;
   for(size_t i = 0; i < splitPoints.size(); ++i)
   {
      if(splitPoints[i] < degree)
      {

         canSwap[splitPoints[i]] = false;
         size_t tmp = _indices[left+splitPoints[i]];
         RCOUT("zostawie" << tmp << " na indeksie "<< splitPoints[i], 11)
         splitPoints[i] = tmp;

      }
      else
      {
         while(!canSwap[index]) index++;
         size_t tmp = _indices[left+splitPoints[i]];
         RCOUT("wstawie" << tmp << " na indeksie "<< left+index,11)
         _indices[left+splitPoints[i]] = _indices[left+index];
         _indices[left+index] = tmp;
         splitPoints[i] = tmp;
         index++;
      }
   }

   for(auto iter = splitPointsRanges_small.begin(); iter != splitPointsRanges_small.end();iter++)
   {
      //Rcout << "rozwazam " << iter->first.i << " i " << iter->first.j <<endl;
      if(find(splitPoints.begin(), splitPoints.end(), iter->first.i) != splitPoints.end() &&
            find(splitPoints.begin(), splitPoints.end(), iter->first.j) != splitPoints.end()) //mamy taki wpis, ze jest to odleglosc miedzy rzeczywiscie wybranymi split pointami
      {
         //Rcout << "wrzucam dla " << iter->first.i << " i " << iter->first.j <<endl;
         splitPointsRanges.insert(*iter);
      }
   }

   return splitPoints;
}

HClustGnatSingleNode* HClustGnatSingle::createNonLeafNode(size_t degree, size_t left, size_t right,const vector<size_t>& splitPoints, const vector<size_t>& boundaries, const vector<size_t>& degrees)
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
      HClustGnatSingleNode* child = buildFromPoints(degrees[i], childLeft, childRight);
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
#if VERBOSE > 11
         Rcout << "szukam/wstawiam dla " <<  splitPoints[mySplitPointIndex] << " i " << splitPoints[j] << endl;
#endif
         auto rangeIterator = splitPointsRanges.find(SortedPoint(splitPoints[mySplitPointIndex], splitPoints[j]));
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
            splitPointsRanges.emplace(SortedPoint(splitPoints[mySplitPointIndex], splitPoints[j]), HClustGnatRange(distances[j], distances[j]));
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
      //Rcout << "cumsum =" <<cumsum << endl;
      boundaries[i] = left+splitPoints.size()+cumsum;
      //Rcout << "boundaries[i] =" << boundaries[i] << endl;
   }
   return boundaries;
}

vector<size_t> HClustGnatSingle::chooseDegrees(size_t left, size_t allPointsCount, const vector<size_t>& boundaries)
{
   /*//left should be after splitpoints!
   vector<size_t> degrees(boundaries.size());
   size_t whole = opts.degree*boundaries.size();
   for(size_t i=0;i<boundaries.size(); ++i)
   {
      size_t elementsCount = i == 0 ? boundaries[i]-left : boundaries[i] - boundaries[i-1];
      degrees[i] = min(max(((elementsCount)/(double)allPointsCount)*whole,(double)opts.minDegree),(double)min(opts.degree*opts.maxTimesDegree, opts.maxDegree));
      Rcout << "degrees[" << i << "]=" << degrees[i] << endl;
   }
   return degrees;*/
   vector<size_t> degrees(boundaries.size());
   for(size_t i=0;i<boundaries.size(); ++i)
   {
      degrees[i] = opts.degree;
   }
   return degrees;
}

HClustGnatSingleNode* HClustGnatSingle::buildFromPoints(size_t degree, size_t left, size_t right)
{
#ifdef GENERATE_STATS
      ++stats.nodeCount;
#endif
   RCOUT("degree = " << degree, 11)
   if(right - left <= opts.candidatesTimes*opts.degree) //@TODO: pomyslec, jaki tak naprawde mamy warunek stopu
   {
      RCOUT("tworze leaf, left= " <<left << " right= " <<right,11)
      return new HClustGnatSingleNode(left, right);
   }
   //Rcout << "split points" << endl;
   vector<size_t> splitPoints = chooseNewSplitPoints(degree, left, right);
#if VERBOSE > 11
   for(size_t i=0;i<splitPoints.size();i++)
   {
      RCOUT("splitpoint[" << i << "]=" << splitPoints[i], 11)
   }
#endif
   //Rcout << "boundaries" << endl;
   printIndices();
   vector<size_t> boundaries = groupPointsToSplitPoints(splitPoints, left, right); //@TODO: dobrze sie zastanowic, gdzie umieszczamy split pointy, aby nie szly w dol, gdzie sa granice!
#if VERBOSE > 11
   for(size_t i=0;i<boundaries.size();i++)
         RCOUT("boundaries[" << i << "]=" << boundaries[i],11);
#endif
   printIndices();
   //@TODO: wybierac degree dziecka, zeby sie roznilo od degree aktualnego, jest w artykule
   vector<size_t> degrees = chooseDegrees(left+degree, right-(left+degree), boundaries);
   //Rcout << "tworze nonleaf" << endl;
   return createNonLeafNode(degree, left, right, splitPoints, boundaries, degrees);
}


void HClustGnatSingle::getNearestNeighborsFromMinRadiusRecursive( HClustGnatSingleNode* node, size_t index,
   size_t clusterIndex, double minR, double& maxR,
   std::priority_queue<HeapNeighborItem>& heap )
{
   // search within (minR, maxR]
   if (node == NULL)
      {
      Rcout << "very strange, node==NULL" << endl;
      return; // this should not happen
      }

#ifdef GENERATE_STATS
   ++stats.nodeVisit;
#endif
   RCOUT("poczatkowe sprawdzanie same cluster",11);
   if (node->sameCluster) {
      if (node->degree == SIZE_MAX) {
         RCOUT("To jest leaf, sprawdzam lewego",11);
         if (ds.find_set(_indices[node->left]) == clusterIndex) return;
      } else {
         RCOUT("To jest non leaf, sprawdzam jego split pointa",11);
         if (node->splitPointIndex != SIZE_MAX && ds.find_set(node->splitPointIndex) == clusterIndex) return; //@TODO: czy na pewno node->splitPointIndex? Jak to wyglada w pierwszym node?
      }
   }
   RCOUT("po sprawdzaniu same cluster z poczatku",11);
   if (node->degree == SIZE_MAX) // leaf
   {
      RCOUT("leaf",11);
      if (node->sameCluster)
      {
         RCOUT("leaf - same cluster",11);
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
            if (heap.size() >= opts.maxNNPrefetch) maxR = heap.top().dist;
         }
      }
      else
      {
         if(node->left == node->right) return; //this is possible and this is nothing wrong. This node is needed, because it has split point. Of course, we could have them in node above...
         RCOUT("leaf - NOT same cluster",11);
         RCOUT("_indices[node->left] = " <<_indices[node->left],11);
         RCOUT("node->left = " << node->left,11);
         RCOUT("node->right = " << node->right,11);
         size_t commonCluster = ds.find_set(_indices[node->left]);
         RCOUT("common cluster sprawdzony",11);
         for (size_t i=node->left; i<node->right; i++)
         {
            size_t currentCluster = ds.find_set(_indices[i]);
            RCOUT("current cluster sprawdzony",11);
            if (currentCluster != commonCluster) commonCluster = SIZE_MAX;
            RCOUT("common cluster != common cluster sprawdzone",11);
            if (currentCluster == clusterIndex) continue;
            RCOUT("common cluster == cluster index sprawdzone",11);
            if (index >= _indices[i]) continue;
            RCOUT("zaraz bede liczyl dystans",11);
            double dist2 = (*_distance)(index, _indices[i]);
            if (dist2 > maxR || dist2 <= minR) continue;
            if (heap.size() >= opts.maxNNPrefetch) {
               if (dist2 < maxR) {
                  while (!heap.empty() && heap.top().dist == maxR) {
                     heap.pop();
                  }
               }
            }
            RCOUT("new object added: " << _indices[i] + 1 << " " << dist2,11)
            heap.push( HeapNeighborItem(_indices[i], dist2) );
            if (heap.size() >= opts.maxNNPrefetch) maxR = heap.top().dist;
         }

         if (commonCluster != SIZE_MAX && ds.find_set(node->splitPointIndex) == commonCluster)
         {node->sameCluster = true;}
      }
      return;
   }
   RCOUT("not leaf",11);
   //Rcout << "I have "<< node->children.size()<<"children" << endl;
   // else // not a leaf
   //1. z artykulu
   vector<bool> shouldVisit(node->degree, true);
   for(size_t i=0;i<node->degree; ++i) //4. z artykulu
   {
      if(shouldVisit[i])
      {
         if(node->children[i]==NULL)
         {
            RCOUT("child is null!",9);
         }
         size_t pi = node->children[i]->splitPointIndex;
         //Rcout << "i got splitPointIndex" << endl;
         //2. z artykulu
         double dist = (*_distance)(node->children[i]->splitPointIndex, index);

         if (ds.find_set(node->children[i]->splitPointIndex) != clusterIndex && index < node->children[i]->splitPointIndex) {
               if (dist <= maxR && dist > minR) {
                  if (heap.size() >= opts.maxNNPrefetch && dist < maxR) {
                     while (!heap.empty() && heap.top().dist == maxR) {
                        heap.pop();
                     }
                  }
                  RCOUT("new object added (splitPoint): " << node->children[i]->splitPointIndex + 1<< " " << dist,11)
                  heap.push( HeapNeighborItem(node->children[i]->splitPointIndex, dist) );
                  if (heap.size() >= opts.maxNNPrefetch) maxR = heap.top().dist;
               }
         }
         // }
         //3. z artykulu
         for(size_t j=0;j<node->degree;++j)
         {
            if(i != j && shouldVisit[j])
            {
               size_t pj = node->children[j]->splitPointIndex;
               //Rcout << "i got splitPointIndex from pj" << endl;
               auto rangeIterator = splitPointsRanges.find(SortedPoint(pi, pj));
               //assert: rangeIterator != splitPointsRanges.end()
               if(rangeIterator == splitPointsRanges.end())
               {
                  RCOUT("distance not found!",11)
               }
               HClustGnatRange range = rangeIterator->second;
               //Rcout << "distance found!" << range.min<< ", " << range.max<<endl;
               double leftRange = dist-maxR;
               double rightRange = dist+maxR;
               if(leftRange <= range.max && range.min <= rightRange) //http://world.std.com/~swmcd/steven/tech/interval.html
               {//they intersect
                  double leftRangeMin  = minR - dist;
                  //double rightRangeMin = minR + dist;
                  if(leftRangeMin >= range.max)
                  {
                     RCOUT("odrzucam ze wzgledu na min R", 8);
                     shouldVisit[j] = false;
                  }
               }
               else
               {//disjoint
                  RCOUT("odrzucam ze wzgledu na max R", 8);
                  shouldVisit[j] = false;
               }
            }
         }

      }
   }

   RCOUT("I have pruned P, i go into my children", 11)

   //5. z artykulu
   for(size_t i=0;i<node->degree; ++i) //4. z artykulu
   {
      if(shouldVisit[i])
      {
         if(node->children[i] == NULL)
         {
            RCOUT("My child is null :(",11);
         }
         getNearestNeighborsFromMinRadiusRecursive(node->children[i], index, clusterIndex, minR, maxR, heap);
      }
   }
   RCOUT("teraz robie rozna magie z same cluster", 11)
   if (node->sameCluster) return;
   RCOUT("sprawdzam, czy dzieci moje sa wszystkie same cluster", 11)
   bool ret = false;
   for(size_t i = 0; i<node->children.size(); ++i)
   {
      if(!node->children[i]->sameCluster)
      {
         ret = true;
         break;
      }
   }
   if(ret) return;
   RCOUT("czy aby wszystkie dzieci moje sa z klastra tego samego?", 11)
   // otherwise check if node->sameCluster flag needs updating
   size_t commonCluster = ds.find_set(node->splitPointIndex);

   for(size_t i = 0; i<node->children.size(); ++i)
   {
      size_t currentCluster = ds.find_set(node->children[i]->splitPointIndex);
      if (currentCluster != commonCluster) return; // not ready yet
   }
   RCOUT("teraz moge spokojnie przypisac samecluster = true", 11)
   node->sameCluster = true;
}

void HClustGnatSingle::FindNeighborTest(size_t index, double R)
{
   std::priority_queue<HeapNeighborItem> heap;
   size_t clusterIndex = ds.find_set(index);
   double _tau = R;
   getNearestNeighborsFromMinRadiusRecursive( _root, index, clusterIndex, minRadiuses[index], _tau, heap );
   Rcout<<"uwaga, koniec, wyniki: " << endl;
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
#if VERBOSE > 11
      Rcout << "minRadius for " << index+1 << " is " <<  minRadiuses[index] << endl;
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
#if VERBOSE > 11
         Rcout << "Nearest for " << i+1 << " is " << hi.index +1 <<" dist=" << hi.dist << endl;
#endif
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
   while(true)
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
         if(i == _n - 1) break;
      }
#if VERBOSE > 7
      if (i % 1024 == 0) Rprintf("\r             %d / %d", i+1, _n);
#endif

      // ASSERT: hhi.index1 < hhi.index2
      HeapNeighborItem hi=getNearestNeighbor(hhi.index1);
      if(hi.index != SIZE_MAX)
      {
#if VERBOSE > 11
         Rcout << "Nearest for " << hhi.index1 +1<< " is " << hi.index+1 <<" dist=" << hi.dist << endl;
#endif
         pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));
      }
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
   //Rcpp::stop("GNAT zbudowany!");
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


