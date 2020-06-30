# fit data

library(tidyverse)
library(cowplot)
source('lib/comp_coef_fns.R')

dir.create('results/figures/coculture')

aij = as.numeric(snakemake@wildcards$aij)
aji = as.numeric(snakemake@wildcards$aji)

outfile = snakemake@output[[1]]

params = read_csv('results/monoculture_parameters.csv')

# initialize dataframe
temps = c('32C','34C','36C','38C','40C')
p_names = c('temp','a_ii','a_jj','a_ij','a_ji','nd','rfd','error')
p = data.frame(matrix(NA, nrow=length(temps),ncol=length(p_names)))
rownames(p) = temps
colnames(p) = p_names

plots = list()

# for each temperature
for(t in temps){

  temp_params = params %>% filter(temp == t)
  attach(temp_params)

  normalized = read_csv(paste0('data/normalized/',temp,'.csv'))
  normalized
  
  # get all data for temperature
  d = normalized %>% gather(Replicate,Growth,EC1:EC3,PP1:PP3,ECco1:ECco3,PPco1:PPco3)
  
  # get EC monoculture data
  ec = filter(d,Replicate %in% c('EC1','EC2','EC3'))
  # get PP monoculture data
  pp = filter(d,Replicate %in% c('PP1','PP2','PP3'))
  # get ECco data
  ecco = filter(d,Replicate %in% c('ECco1','ECco2','ECco3'))
  # get PPco coculture data
  ppco = filter(d,Replicate %in% c('PPco1','PPco2','PPco3'))
  # get all co-culture data
  dat = filter(d,Replicate %in% c('ECco1','ECco2','ECco3','PPco1','PPco2','PPco3'))
  
  a_ii_start = a_jj_start = a_ij_start = a_ji_start = 0.5
  
  if(temp == '32C'){
    dat = dat[dat$Time > 6,]
  }
  if(temp == '34C'){
    dat = dat[dat$Time > 3.5,]
  }
  if(temp == '36C'){
    dat = dat[dat$Time > 5,]
  }
  if(temp == '38C'){
    dat = dat[dat$Time > 5,]
  }
  if(temp == '40C'){
    dat = dat[dat$Time > 5,]
  }

  fit <- try(optim(par=c(aij,aji,ecco_n0,ppco_n0,a_ii_start,a_jj_start),
                   fn=ofun,control = list(maxit = 1000),
                   lower=c(aij,aji,0,0,0,0), upper=c(aij+10e-5,aji+10e-5,Inf,Inf,Inf,Inf),method='L-BFGS-B'))
  if(inherits(fit,"try-error")){ 
     p[temp,] = c(temp, NA, NA, NA, NA, NA, NA, NA)
     next 
  }

  fitdat <- data.frame(Time=dat$Time,
                       curve=comp.coeff.curve(fit$par[1],fit$par[2],fit$par[3],fit$par[4],fit$par[5],fit$par[6])[,1:2],
                       dat=dat$Growth,Replicate=dat$Replicate)
  
  plots[[temp]] = ggplot(fitdat,aes(x=Time,y=dat)) +
    geom_point(aes(col=Replicate))+geom_line(aes(y=curve.E)) +
    geom_line(aes(y=curve.P)) + ggtitle(temp) + theme_bw() + ylab('Growth (r.f.u.)') +
    theme(legend.position="none",plot.margin = margin(6, 0, 6, 0))
  
  # fitted
  a_ij = fit$par[1]
  a_ji = fit$par[2]
  a_ii = fit$par[5]
  a_jj = fit$par[6]
  
  # calculate ND
  nd = 1 - sqrt((a_ij*a_ji)/(a_ii*a_jj))
  
  # calculate RFD
  rfd = sqrt((a_ij*a_ii)/(a_ji*a_jj))
  
  # error
  err = fit$value
  
  # save all to dataframe
  p[temp,] = c(temp, a_ii, a_jj, a_ij, a_ji, nd, rfd, err)
  
}

write_csv(p,outfile)

if(length(plots) != 0){
  plot_fitted(plots,aij,aji)
}
