library(tidyverse)
library(magrittr)

BINWIDTH <- 0.05

neurdf <- read_csv("./trialspikedata.csv", guess_max = 100000)
trialdf<- neurdf %>% filter(trial>0, spiketime <= completiontime, spiketime >= trialstart)

make_counts <- function(times, start, end, width){
  cuts <- seq(start, end, width)
  if(cuts[length(cuts)] < end){
    cuts <- c(cuts,cuts[length(cuts)]+width)
  }  
  counts <- rep(0, length(cuts)-1)
  b <- .bincode(times, cuts, right=FALSE)
  ub <- unique(b)
  counts[ub] <- tabulate(b)[ub]
  return(counts)
}


trialdata <- trialdf %>% select(-spiketime) %>% distinct()

binned_df <- trialdf %>% group_by(region, neuron_number, trial) %>% 
  group_modify(~ { data_frame(bincount =make_counts(.x$spiketime, .x$trialstart[1], .x$completiontime[nrow(.x)], BINWIDTH),
                              binnumb=1:length(bincount), binstart=(binnumb-1)*BINWIDTH + .x$trialstart[1], binend=(binnumb)*BINWIDTH + .x$trialstart[1])})


binned_df %<>% left_join(trialdata, by=c("region", "neuron_number", "trial"))
