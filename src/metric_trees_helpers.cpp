#include "metric_trees_helpers.h"

using namespace Rcpp;
using namespace std;

ostream& operator<< (ostream& os, const Point& obj) {
   os << obj.i<< " " << obj.j;
   return os;
}

istream& operator>> (istream& is, Point& obj) {
   is >> obj.i;
   is >> obj.j;
   return is;
}

vector<RObject> createStdVectorOfRobjects(List listobj)
{
   int n = listobj.size();
   vector<RObject> points(n);
   for(int i=0;i<n;i++)
   {
      RObject obj = listobj[i];
      R_PreserveObject(obj);
      points[i] = obj;
   }
   return points;
}
