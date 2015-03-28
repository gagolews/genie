print.Vector <- function(vec)
{
   print(as.list.Vector(vec))
}

format.Vector <- function(vec, ...)
{
   format(as.list.Vector(vec), ...)
}

# "[.Vector" <- function(vec, i)
# {
#    vector_at(vec,i)
# }
# 
# "[.Vector" <- function(vec, i, value)
# {
#    vector_at(vec,i) = value
# }
# 
"vector_at<-" <- function(vec, i, value)
{
   vector_set_at(vec,i,value)
   vec
}