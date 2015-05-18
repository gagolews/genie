#' @title Queue
#'
#' @description
#' Queue is an abstract data type which allows for pushing elements
#' and popping them in normal (first in-first out) order.
#'
#' @rdname queue
#' @name Queue
invisible(NULL)

#' @rdname queue
#' @details
#' \code{print.Queue} prints the contents of the queue
#' (coerced to a list object) on the console.
#'
#' @return
#' \code{print.Queue} returns the queue coerced to a list, invisibly.
#'
#' @param ... other arguments passed to or from other methods
print.Queue <- function(queue, ...)
{
   print(as.list.Queue(queue), ...)
}

#' @rdname queue
#' @details
#' \code{format.Queue} pretty-prints the contents of the queue
#' (coerced to a list object) on the console.
#'
#' @return
#' \code{format.Queue} returns the character representation
#' of the objects in the queue
#'
#' @param ... other arguments passed to or from other methods
format.Queue <- function(queue, ...)
{
   format(as.list.Queue(queue), ...)
}

#' @rdname queue
#' @details
#' \code{str.Queue} compactly displays the contents of the queue
#' (coerced to a list object) on the console.
#'
#' @return
#' \code{str.Queue} does not return anything interesting.
str.Queue <- function(queue, ...)
{
   str(as.list.Queue(queue), ...)
}