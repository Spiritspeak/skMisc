

#' @name RateLimit
#' @title Function call rate limiting
#'
#' @param period Time period in seconds within which the rate of function calls must be limited.
#' @param times Maximum number of times this queue can be utilized within the designated \code{period}.
#' @param padding Additional amount of time in seconds to wait until the script continues, 
#' when the rate limit has been exceeded.
#' @param queue The rate limit object.
#'
#' @returns
#' @export
#'
#' @examples
#' h<-defineRateLimit(period=60,times=30)
#' for(i in 0:30){
#'   message(i)
#'   awaitRateLimit(h)
#' }
#' 
defineRateLimit <- function(period, times, padding=0){
  env<-new.env()
  env$period <- period
  env$times <- times
  env$padding <- padding
  env$events <- c()
  structure(env, class="RateLimitQueue")
}


#' @export
#' @rdname RateLimit
awaitRateLimit <- function(queue){
  tds <- as.numeric(difftime(Sys.time(), queue$events, units="secs"))
  queue$events <- queue$events[tds < queue$period]
  if(length(queue$events) >= queue$times){
    Sys.sleep(max(0, queue$period - max(tds) + queue$padding))
  }
  tds2 <- as.numeric(difftime(Sys.time(), queue$events, units="secs"))
  queue$events <- c(Sys.time(), queue$events[tds2 < queue$period])
}

