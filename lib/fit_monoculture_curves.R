# fit monoculture curves

library(tidyverse)
source('lib/comp_coef_fns.R')

temps = c('32C','34C','36C','38C','40C')
# p_names = c('temp','ec_r','pp_r','ec_k','pp_k')
# rownames(p) = temp
# colnames(p) = p_names

p = data.frame(t(sapply(temps, function(temp){
  # path to normalized data
  normalized = read_csv(paste0('data/normalized/',temp,'.csv'))
  
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

  # get EC parameters
  ec_k = est_k(ec$Time,ec$Growth) #ec$Growth[nrow(ec)]
  ec_r = est_r(ec$Time,ec$Growth)
  # get PP parameters
  pp_k = est_k(pp$Time,pp$Growth) #pp$Growth[nrow(pp)]
  pp_r = est_r(pp$Time,pp$Growth)
  # get ECco parameters
  ecco_k = est_k(ecco$Time,ecco$Growth)
  ecco_n0 = est_n0(ecco$Time,ecco$Growth)
  # get PPco parameters
  ppco_k = est_k(ppco$Time,ppco$Growth)
  ppco_n0 = est_n0(ppco$Time,ppco$Growth)
  return(c(temp=temp,ec_k=ec_k,ec_r=ec_r,pp_k=pp_k,pp_r=pp_r,ecco_k=ecco_k,ecco_n0=ecco_n0,ppco_k=ppco_k,ppco_n0=ppco_n0))
})))

write_csv(p,'results/monoculture_parameters.csv')
