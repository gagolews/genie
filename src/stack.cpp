#include <Rcpp.h>
using namespace Rcpp;
#include <deque>

class S : public std::deque<SEXP> {
   public: ~S() {
//      Rcout << "DESCRUCT!\n";
      while (!this->empty()) {
         SEXP obj = this->front();
         this->pop_front();
         R_ReleaseObject(obj);
      }
   }

   static const char* ClassName;
};


const char* S::ClassName = "Stack";


void checkIsStackClass(XPtr< S >& _stack)
{
   if (strcmp(_stack.attr("class"), S::ClassName))
      stop("not a Stack object");
}


//' @rdname stack
//' @details
//' \code{stack_create} creates a new stack.
//'
//' @return
//' \code{stack_create} returns a new, empty stack.
// [[Rcpp::export]]
SEXP stack_create() {
   S* stack = new S();
   XPtr< S > retval =  XPtr< S >(stack, true);
   retval.attr("class") = S::ClassName;
   return retval;
}


//' @rdname stack
//' @details
//' \code{as.list.Stack} converts a given stack object to an R list.
//'
//' @return
//' \code{as.list.Stack} returns an R list object.
//' @param stack a stack object
// [[Rcpp::export("as.list.Stack")]]
List stack_as_list(SEXP stack) {
   XPtr< S > _stack = Rcpp::as< XPtr< S > > (stack);
   checkIsStackClass(_stack);
   int n = (*_stack).size();
   List retval(n);
   std::deque<SEXP>::iterator it = (*_stack).begin();
   for (int i=0; i < n; ++i)
      retval[i] = *(it++);
   return retval;
}

//' @rdname stack
//' @details
//' \code{stack_empty} determines if a given stack is empty or not.
//'
//' @return
//' \code{stack_empty} returns a single logical value.
// [[Rcpp::export]]
bool stack_empty(SEXP stack) {
   XPtr< S > _stack = Rcpp::as< XPtr< S > > (stack);
   checkIsStackClass(_stack);
   return (*_stack).empty();
}


//' @rdname stack
//' @details
//' \code{stack_push} pushes a given object to the top of the stack.
//'
//' @return
//' \code{stack_push} does not return anything interesting.
//'
//' @param obj an R object
// [[Rcpp::export]]
void stack_push(SEXP stack, SEXP obj) {
   XPtr< S > _stack = Rcpp::as< XPtr< S > > (stack);
   checkIsStackClass(_stack);
   R_PreserveObject(obj);
   (*_stack).push_front(obj);
}


//' @rdname stack
//' @details
//' \code{stack_pop} pops an object from the top of the stack.
//'
//' @return
//' \code{stack_pop} returns an object at the top of te stack or
//' throws an error if the stack is empty.
// [[Rcpp::export]]
SEXP stack_pop(SEXP stack) {
   XPtr< S > _stack = Rcpp::as< XPtr< S > > (stack);
   checkIsStackClass(_stack);
   if ((*_stack).empty())
      stop("empty stack");
   SEXP obj = (*_stack).front();
   (*_stack).pop_front();
   R_ReleaseObject(obj);
   return obj;
}
