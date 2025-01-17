

defineRateLimit <- function(period, times, padding=0){
  env<-new.env()
  env$period <- period
  env$times <- times
  env$padding <- padding
  env$events <- c()
  structure(env, class="RateLimitQueue")
}

#h<-defineRateLimit(period=60,times=30)
#awaitRateLimit(h)

awaitRateLimit <- function(queue){
  tds <- as.numeric(difftime(Sys.time(), queue$events, units="secs"))
  queue$events <- queue$events[tds < queue$period]
  if(length(queue$events) >= queue$times){
    Sys.sleep(max(0, queue$period - max(tds) + queue$padding))
  }
  tds2 <- as.numeric(difftime(Sys.time(), queue$events, units="secs"))
  queue$events <- c(Sys.time(), queue$events[tds2 < queue$period])
}


