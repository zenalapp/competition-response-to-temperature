# functions to find competition coefficients for cocultures using the lotka voltera equation

# load packages ----
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(deSolve))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(growthcurver))

# estimate parameters ----
# estimate k
est_k = function(time,data){
  ests = SummarizeGrowth(time,data)
  k = ests$vals$k
  return(k)
}

# estimate r
est_r = function(time,data){
  ests = SummarizeGrowth(time,data)
  #plot(ests)
  r = ests$vals$r
  return(r)
}

# estimate n0
est_n0 = function(time,data){
  ests = SummarizeGrowth(time,data)
  n0 = ests$vals$n0
  return(n0)
}

# model functions ----
# model (lotka voltera equations)
comp.coeff.model <- function(t,x,params) {
  #extract the state variables
  E <- x[1]
  P <- x[2]
  ec_r = ec_r
  ec_k = ec_k
  pp_r = pp_r
  pp_k = pp_k
  #extract the parameters
  a_ij <- params["a_ij"]
  a_ji <- params["a_ji"]
  a_ii = params["a_ii"]#1/ec_k
  a_jj = params["a_jj"]#1/pp_k 
  #now code the model equations
  # dEdt <- ec_r*E*(1-a_ii*E/ec_k-a_ij*P/ec_k)
  # dPdt <- pp_r*P*(1-a_jj*P/pp_k-a_ji*E/pp_k)
  dEdt <- ec_r*E*(1-(a_ii*E)-(a_ij*P))
  dPdt <- pp_r*P*(1-(a_jj*P)-(a_ji*E))
  #combine results into a single vector
  dxdt <- c(dEdt,dPdt)
  #return results as a list!
  list(dxdt)
}

# get curve fitted to model
comp.coeff.curve <- function(a_ij,a_ji,E0,P0,a_ii,a_jj) {
  xstart <- c(E=E0, P=P0)
  params <- c(a_ij=a_ij,a_ji=a_ji,E0=E0,P0=P0,a_ii=a_ii,a_jj=a_jj)
  times <- unique(dat$Time)
  traj <- ode(func = comp.coeff.model,y=xstart,
              times = times,parms = params,
              maxsteps = 10000, method = "rk4")
  y <- traj[,-1,drop=FALSE]
  return(y)
}

# error function
error <- function(a_ij,a_ji,E0,P0,a_ii,a_jj) {
  Y <- comp.coeff.curve(a_ij,a_ji,E0,P0,a_ii,a_jj)
  ecdat = filter(dat,Replicate %in% c('ECco1','ECco2','ECco3'))
  ppdat = filter(dat,Replicate %in% c('PPco1','PPco2','PPco3'))
  sum((cbind(((ecdat[,'Growth']-rep(Y[,1],3))/ecdat[,'Growth'])^2,
             ((ppdat[,'Growth']-rep(Y[,2],3))/ppdat[,'Growth'])^2)))
}

# calculate error
ofun <- function(par) {
  error(a_ij=par[1],a_ji=par[2],E0=par[3],P0=par[4],a_ii=par[5],a_jj=par[6])
}

# fit function ----
plot_fitted = function(plots,aij,aji){
  pgrid <- plot_grid(plots,
    align = 'vh'
  )
  # legend
  legend = get_legend(
    # create some space to the left of the legend
    plots[[1]] + theme(legend.box.margin = margin(0, 0, 0, 0))
  )
  p <- plot_grid(pgrid, legend, 
                  rel_widths = c(3, .5), nrow=1)
  print(p)
  ggsave(filename=paste0('results/figures/coculture/fit_coculture_aij_',aij,'_aji_',aji,'.pdf'),plot=p)
}

#---------------------------------------------
# OLD FUNCTIONS

# plot original data
plot_all_data = function(mat,time='Time',species1='ECco',species2='PPco',main=NULL){
  ggplot(data=mat,mapping=aes_string(x=time,y=species1,col=shQuote(species1)))+geom_point()+geom_point(aes_string(y=species2,col=shQuote(species2))) + ylab('Growth (OD600)') + ggtitle(main) + theme_classic()
}

# estimate k and r
est_k_r = function(time,data){
  ests = SummarizeGrowth(time,data)
  #plot(ests)
  k = ests$vals$k
  r = ests$vals$r
  n0 = ests$vals$n0
  return(c(k=k,r=r,n0=n0))
}

# # plot fitted
# plot_fitted = function(dat,fit,main=NULL){
#   Y <- comp.coeff.curve(fit$par[1],fit$par[2],fit$par[3],fit$par[4])
#   fitdat <- as.data.frame(Y)
#   #fitdat <- as.data.frame(Y[seq(1,nrow(Y),2),])
#   dat <- mutate(dat,fitdatE=fitdat$E,fitdatP=fitdat$P)
#   ggplot(data=dat,mapping=aes(x=Time,y=ECco))+
#     geom_point(aes(col='EC'))+geom_line(aes(y=fitdatE))+
#     geom_point(aes(y=PPco,col='PP'))+
#     geom_line(aes(y=fitdatP)) + ylab('Growth (OD600)') + ggtitle(main) + theme_bw()
# }

# get breakpoint using broken stick
get_breakpoint = function(dat){
  linear_ec = lm(ECco ~ Time, dat)
  linear_pp = lm(PPco ~ Time, dat)
  # broken stick
  bp = median(sapply(1:10, function(x){
    bs_ec = segmented(linear_ec,psi = c(1,10))
    bp_ec = bs_ec$psi[1,2]
    bs_pp = segmented(linear_pp,psi = c(1,10))
    bp_pp = bs_pp$psi[1,2]
    bp = min(bp_ec,bp_pp)
    slope_diff_ec = bs_ec$coefficients[3] - bs_ec$coefficients[2]
    slope_diff_pp = bs_pp$coefficients[3] - bs_pp$coefficients[2]
    if(slope_diff_ec < -1.5 | slope_diff_pp < -1.5) bp = 0 # one lag phase starts from beginning (ideally)
    return(bp)
  }))
  return(bp)
}


# get breakpoint using broken stick
get_exp_phase = function(time,growth){
  # get fitted curve
  fit = SummarizeGrowth(time,growth)
  fitted = fit$model$m$fitted()
  # make dataframe
  dat = data.frame(time=time,fitted=fitted)
  # linear model to start
  linear = lm(fitted ~ time, dat)
  # broken stick
  bs = segmented(linear, psi=c(1,10))
  first_break = bs$psi[1,2]
  second_break = bs$psi[2,2]
  #plot(fit)
  #abline(v=second_break+2)
  return(c('start_exp'=first_break,'end_exp'=second_break))
}


# dealing with non-identifiability ----

get_error = function(aij,aji){
  fit <- try(optim(par=c(aij,aji,ecco_n0,ppco_n0,a_ii_start,a_jj_start),
                   fn=ofun,control = list(maxit = 1000),
                   lower=c(aij,aji,0,0,0,0), upper=c(aij+10e10,aji+10e10,Inf,Inf,Inf,Inf),method='L-BFGS-B'))
  if(inherits(fit,"try-error")) return(NA)
  return(fit$value)
}

