#' @title Vector
#'
#' @description
#' Vectors are sequence containers representing arrays that can change in size.
#'
#' @rdname vector
#' @name Vector
invisible(NULL)

#' @rdname vector
#' @details
#' \code{print.Vector} prints the contents of the vector
#' (coerced to a list object) on the console.
#'
#' @return
#' \code{print.Vector} returns the vector coerced to a list, invisibly.
#'
#' @param ... other arguments passed to or from other methods
print.Vector <- function(vec, ...)
{
   print(as.list.Vector(vec), ...)
}

#' @rdname vector
#' @details
#' \code{format.Vector} pretty-prints the contents of the vector
#' (coerced to a list object) on the console.
#'
#' @return
#' \code{format.Vector} returns the character representation
#' of the objects in the vector.
#'
#' @param ... other arguments passed to or from other methods
format.Vector <- function(vec, ...)
{
   format(as.list.Vector(vec), ...)
}

#' @rdname vector
#' @details
#' \code{[.Vector} returns an object at position i in the vector.
#'
#' @return
#' \code{[.Vector} returns a single R object.
#'
#' @param i position i in the vector
`[.Vector` <- function(vec, i)
{
   vector_at(vec,i)
}

#' @rdname vector
#' @details
#' \code{[<-.Vector} sets an object value at position i in the vector.
#'
#' @return
#' \code{[<-.Vector} returns a modified vector
#'
#' @param i position i in the vector
#' @param value object to set
`[<-.Vector` <- function(vec, i, value)
{
   vector_set_at(vec,i,value)
   vec
}

#' @rdname vector
#' @details
#' \code{vector_at<-} sets an object value at position i in the vector.
#'
#' @return
#' \code{vector_at<-} returns a modified vector
#'
#' @param i position i in the vector
#' @param value object to set
"vector_at<-" <- function(vec, i, value)
{
   vector_set_at(vec,i,value)
   vec
}
