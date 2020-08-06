library(tidyverse)

theme_set(theme_bw())

all_results <- read_csv('results/all_results.csv')

p <- all_results %>% 
  filter(!is.na(nd) & !is.na(rfd) & error < 100 & !is.infinite(nd) & !is.infinite(rfd) & rfd < 50) %>% 
  ggplot(aes(x=nd,y=rfd,col=temp,alpha=1/log10(error))) + geom_jitter() + facet_wrap(vars(temp))

ggsave(p, file = 'results/figures/nd_rfd/nd_rfd.pdf')
