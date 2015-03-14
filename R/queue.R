print.Queue <- function(queue)
{
   print(as.list.Queue(queue))
}

format.Queue <- function(queue, ...)
{
   format(as.list.Queue(queue), ...)
}
