# plot error

library(tidyverse)

dat <- read_csv('results/all_results.csv')

p <- dat %>% select(temp,a_ij,a_ji,error) %>% 
    ggplot(aes(x=a_ij,y=a_ji,z=log10(error))) +
    #geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
    geom_point(aes(color=log10(error))) + scale_color_gradient(low='cornflowerblue',high='white') + 
    geom_hline(yintercept=0) + geom_vline(xintercept=0) +
    facet_wrap(~temp) + theme_bw()  #+
    #xlim(c(-1,1)) + ylim(c(-1,1))

hist <- dat %>% select(temp,a_ij,a_ji,error) %>% filter(a_ij >= 0 & a_ji >= 0)  %>% 
    ggplot(aes(x=error)) +
    geom_vline(xintercept=100,col='red') + #geom_vline(xintercept=0) + 
    geom_histogram(aes(color=log(error)),bins=100) + scale_color_gradient(low='cornflowerblue',high='white') + 
    facet_wrap(~temp) + theme_bw()  #+
    #xlim(c(-1,1)) + ylim(c(-1,1))

p2 <- dat %>% select(temp,a_ij,a_ji,error) %>% filter(a_ij >= 0 & a_ji >= 0)  %>%
    ggplot(aes(x=a_ij,y=a_ji,z=log10(error))) +
    #geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
    geom_point(aes(color=log10(error))) + scale_color_gradient(low='cornflowerblue',high='white') + 
    #geom_hline(yintercept=0) + geom_vline(xintercept=0) +
    facet_wrap(~temp) + theme_bw() 

ggsave('results/figures/error/error.pdf', p)
ggsave('results/figures/error/error_0-5.pdf', p2)
ggsave('results/figures/error/error_hist_0-5.pdf',hist)
