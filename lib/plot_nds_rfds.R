library(tidyverse)

theme_set(theme_bw())

all_results <- read_csv('results/all_results.csv')

p <- all_results %>% 
  filter(!is.na(nd) & !is.na(rfd) & 
           #error < 100 & 
           !is.infinite(nd) & !is.infinite(rfd) & 
           rfd < 10 & a_ij > 0 & a_ji > 0) %>% 
  ggplot(aes(x=nd,y=rfd,col=temp,alpha=1/log10(error))) + geom_jitter() #+ facet_wrap(vars(temp))

ggsave(p, file = 'results/figures/nd_rfd/nd_rfd.pdf')

# all_results %>% 
#   filter(!is.na(nd) & !is.na(rfd) & 
#            error < 50 & 
#            !is.infinite(nd) & !is.infinite(rfd) & 
#            rfd < 10 & a_ij > 0 & a_ji > 0) %>% 
#   group_by(temp) %>% summarize(nd_mean=mean(nd),rfd_mean=mean(rfd)) %>% ggplot(aes(x=nd_mean,y=rfd_mean,col=temp)) + geom_point(size=5)

#all_results %>% group_by(temp) %>% filter(error == min(error,na.rm=T))