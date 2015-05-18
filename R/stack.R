#' @title Stack
#'
#' @description
#' Stack is an abstract data type which allows for pushing elements
#' and popping them in reverse (last in-first out) order.
#'
#' @rdname stack
#' @name Stack
invisible(NULL)


#' @rdname stack
#' @details
#' \code{print.Stack} prints the contents of the stack
#' (coerced to a list object) on the console.
#'
#' @return
#' \code{print.Stack} returns the stack coerced to a list, invisibly.
#'
#' @param ... other arguments passed to or from other methods
print.Stack <- function(stack, ...)
{
   print(as.list.Stack(stack), ...)
}

#' @rdname stack
#' @details
#' \code{format.Stack} pretty-prints the contents of the stack
#' (coerced to a list object) on the console.
#'
#' @return
#' \code{format.Stack} returns the character representation
#' of the objects in the stack.
#'
#' @param ... other arguments passed to or from other methods
format.Stack <- function(stack, ...)
{
   format(as.list.Stack(stack), ...)
}


#' @rdname stack
#' @details
#' \code{str.Stack} compactly displays the contents of the stack
#' (coerced to a list object) on the console.
#'
#' @return
#' \code{str.Stack} does not return anything interesting.
str.Stack <- function(stack, ...)
{
   str(as.list.Stack(stack), ...)
}
