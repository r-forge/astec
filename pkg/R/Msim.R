Msim <- function(parameters, stats_output) {
  .C("metacommunity_model", parameters = c(1,0.5,1,0.5), stats_output=array(0,4))
}

