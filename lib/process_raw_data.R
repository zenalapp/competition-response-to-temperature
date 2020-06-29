# process raw data
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(readxl))

#### Reading & plotting raw data ####

get_data = function(df,inames,start,end){
  df = df[start:end,]
  names(df) = unlist(df[inames,])
  df = df[,!apply(df,2,function(x) all(x==''))]
  # remove water wells
  df = df[,-c(3,13,14)]
  # change colnames
  names(df) = c('Time','Temp','EC1','EC2','EC3','ECPP1','ECPP2','ECPP3','PP1','PP2','PP3')
  df$Time = period_to_seconds(hms(df$Time))/60/60
  df[] <- sapply(df, function(x) as.numeric(as.character(x)))
  df = df[df$Time != 0,]
  return(df)
}

get_data_36 = function(df,inames,start,end){
  # water: 4,14,15,16,26,27,28,38,39
  df = df[start:end,2:ncol(df)]
  names(df) = unlist(df[inames,])
  df = df[,!apply(df,2,function(x) all(x==''))]
  # keep only certain wells
  df = df[,c(1,2,4,5,6,21,22,23,28,29,30)]
  # change colnames
  names(df) = c('Time','Temp','EC1','EC2','EC3','ECPP1','ECPP2','ECPP3','PP1','PP2','PP3')
  df$Time = period_to_seconds(hms(df$Time))/60/60
  df[] <- sapply(df, function(x) as.numeric(as.character(x)))
  df = df[df$Time != 0,]
  return(df)
}

read_raw_data = function(dir = 'data/raw/'){
  # file paths
  raw_files = list.files(dir,full.names = T)
  # initialize lists for data
  ods = list()
  ecf = list()
  ppf = list()
  # loop over each file
  for(i in raw_files){
    # temperature
    temp = gsub('.csv|.xlsx','',strsplit(i,' ')[[1]][3])
    if(temp == '36C'){
      dat = read.csv(i)
      ods[[temp]] = get_data_36(dat,50,51,153)
      ecf[[temp]] = get_data_36(dat,199,200,302)
      ppf[[temp]] = get_data_36(dat,348,349,451)
    }else{
      # read in data
      dat = read.csv(i)
      # ods
      ods[[temp]] = get_data(dat,50,51,195)
      # ec fluorescence (510_540)
      ecf[[temp]] = get_data(dat,199,200,344)
      # pp fluorescence (585_615)
      ppf[[temp]] = get_data(dat,348,349,493)
    }
  }
  return(list(ods=ods,ecf=ecf,ppf=ppf))
}

plot_raw = function(dat,dat_type){
  temps = names(dat)
  plots = list()
  for(i in 1:length(dat)){
    temp = temps[i]
    vals = dat[[i]]
    if(dat_type == 'Normalized'){
      plots[[temp]] = vals %>%
        gather(Replicate,Growth,EC1:EC3,PP1:PP3,ECco1:ECco3,PPco1:PPco3) 
    }else{
      plots[[temp]] = vals %>% gather(Replicate,Growth,EC1:EC3,PP1:PP3,ECPP1:ECPP3) 
    }
    plots[[temp]] = plots[[temp]] %>% ggplot(aes(x=Time,y=Growth,col=Replicate)) + 
            geom_point() + theme_bw() + ggtitle(dat_type)
  }
  #plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]])
  return(plots)
}

plot_plots = function(dat,norm_dat){
  plots_ods = plot_raw(dat$ods,'ODs')
  plots_ecf = plot_raw(dat$ecf,'E. coli fluorescence')
  plots_ppf = plot_raw(dat$ppf,'P. putida fluorescence')
  plots_norm = plot_raw(norm_dat,'Normalized')
  for(i in names(plots_ods)){
    # title
    title = ggdraw() + 
      draw_label(i, fontface = 'bold', x = 0, hjust = 0) +
      theme(plot.margin = margin(0, 0, 0, 0)
    )
    # plots
    plots <- plot_grid(
      plots_ecf[[i]] + theme(legend.position="none",plot.margin = margin(6, 0, 6, 0)),
      plots_ppf[[i]] + theme(legend.position="none",plot.margin = margin(6, 0, 6, 0)),
      plots_ods[[i]] + theme(legend.position="none",plot.margin = margin(6, 0, 6, 0)),
      plots_norm[[i]] + theme(legend.position="none",plot.margin = margin(6, 0, 6, 0)),
      align = 'vh'
    )
    # legend
    legend = get_legend(
      # create some space to the left of the legend
      plots_ecf[[i]] + theme(legend.box.margin = margin(0, 0, 0, 0))
    )
    # legend2
    legend2 = get_legend(
      # create some space to the left of the legend
      plots_norm[[i]] + theme(legend.box.margin = margin(0, 0, 0, 0)) + 
        guides(color=guide_legend(title="Normalized\nReplicate"))
    )
    # combine title and plots
    comb = plot_grid(title, plots, rel_widths = c(0.2,3), rel_heights = c(0.1,1), ncol=1)
    # plot
    print(plot_grid(comb, legend, legend2, 
              rel_widths = c(3, .5, .5), nrow=1))
  }
}

#### Processing data ####

flatten_stationary = function(dat,trunc_flatten){
  ecf_dat = dat$ecf
  ppf_dat = dat$ppf
  temps = names(ecf_dat)
  for(temp in temps){
    ecf = ecf_dat[[temp]]
    ppf = ppf_dat[[temp]]
    # flatten out stationary phase
    ec_reps = c('EC1','EC2','EC3')
    for(j in ec_reps){
        max_fl = trunc_flatten$fl_flat[trunc_flatten$Temp == temp & trunc_flatten$Replicate == j]
        if(is.na(max_fl)) next
        end_exp = min(ecf$Time[ecf[,j] > max_fl])
        ecf[,j][ecf$Time>=end_exp] = max_fl
    }
      pp_reps = c('PP1','PP2','PP3')
      for(j in pp_reps){
        max_fl = trunc_flatten$fl_flat[trunc_flatten$Temp == temp & trunc_flatten$Replicate == j]
        if(is.na(max_fl)) next
        end_exp = min(ppf$Time[ppf[,j] > max_fl])
        ppf[,j][ppf$Time>=end_exp] = max_fl
      }
      co_reps = c('ECPP1','ECPP2','ECPP3')
      for(j in co_reps){
        r = gsub('ECPP','',j)
        ecco = paste0('ECco',r)
        ppco = paste0('PPco',r)
        max_fl_ec = trunc_flatten$fl_flat[trunc_flatten$Temp == temp & trunc_flatten$Replicate == ecco]
        if(!is.na(max_fl_ec)){
          end_exp_ec = min(ecf$Time[ecf[,j] > max_fl_ec])
        }else{
          max_fl_ec = ecf[nrow(ecf),j]
          end_exp_ec = ecf$Time[nrow(ecf)]
        }
        max_fl_pp = trunc_flatten$fl_flat[trunc_flatten$Temp == temp & trunc_flatten$Replicate == ppco]
        if(!is.na(max_fl_pp)){
          end_exp_pp = min(ppf$Time[ppf[,j] > max_fl_pp])
        }else{
          max_fl_pp = ppf[nrow(ppf),j]
          end_exp_pp = ppf$Time[nrow(ppf)]
        }
        end_exp = min(end_exp_ec,end_exp_pp)
        ecf[,j][ecf$Time>=end_exp] = ecf[,j][ecf$Time==ecf$Time[ecf$Time >= end_exp][1]]
        ppf[,j][ppf$Time>=end_exp] = ppf[,j][ppf$Time==ppf$Time[ppf$Time >= end_exp][1]]
      }
      ecf_dat[[temp]] = ecf
      ppf_dat[[temp]] = ppf
  }
  return(list(ods=dat$ods,ecf=ecf_dat,ppf=ppf_dat))
}

normalize_fluorescence = function(dat){
  all_data = list()
  temps = names(dat$ods)
  for(temp in temps){
    ods = dat$ods[[temp]]
    ecf = dat$ecf[[temp]]
    ppf = dat$ppf[[temp]]
    # initialize dataframe
    all_in_ecf = data.frame(matrix(NA,nrow=nrow(ods),ncol=ncol(ods)+3))
    colnames(all_in_ecf) = c('Time','Temp','EC1','EC2','EC3','ECco1','ECco2','ECco3','PPco1','PPco2','PPco3','PP1','PP2','PP3')
    # add time and temp
    all_in_ecf[,c('Time','Temp')] = ecf[,c('Time','Temp')]
    # convert PP fluorescence units to EC fluorescence units (monoculture)
    # get max OD for EC monoculture
    max_ecod_mono = mean(unlist(ods[nrow(ods),c('EC1','EC2','EC3')]))
    # get max OD for PP monoculture
    max_ppod_mono = mean(unlist(ods[nrow(ods),c('PP1','PP2','PP3')]))
    # get max fluorescence for EC monoculture
    max_ecf_mono = mean(unlist(ecf[nrow(ecf),c('EC1','EC2','EC3')]))
    # get max fluorescence for PP monoculture
    max_ppf_mono = mean(unlist(ppf[nrow(ppf),c('PP1','PP2','PP3')]))
    if(temp == '32C'){
      max_ecod_mono_32 = max_ecod_mono
      max_ecf_mono_32 = max_ecf_mono
      max_ppod_mono_32 = max_ppod_mono
    }
    # get normalization value (all relative to 32C)
    norm_ec = max_ecf_mono_32*max_ecod_mono/max_ecod_mono_32/max_ecf_mono
    norm_pp = max_ecf_mono_32*max_ppod_mono/max_ecod_mono_32/max_ppf_mono
    
    # add EC monoculture and coculture data to dataframe
    all_in_ecf[,c('EC1','EC2','EC3','ECco1','ECco2','ECco3')] = 
      ecf[,c('EC1','EC2','EC3','ECPP1','ECPP2','ECPP3')]*norm_ec/10000
    # add PP monoculture and coculture data to dataframe
    all_in_ecf[,c('PPco1','PPco2','PPco3','PP1','PP2','PP3')] = 
      ppf[,c('ECPP1','ECPP2','ECPP3','PP1','PP2','PP3')]*norm_pp/10000
    # tidy data
    #all_in_ecf = all_in_ecf %>% gather(Replicate,Growth,EC1:EC3,PP1:PP3,ECco1:ECco3,PPco1:PPco3)
    #ods = ods %>% gather(Replicate,Growth,EC1:EC3,PP1:PP3,ECPP1:ECPP3)
    all_data[[temp]] = all_in_ecf
  }
  return(all_data)
}

# 
# ###################
# 
# process_raw_data = function(){
#   
#   trunc_flatten = read.csv('trunc_flatten.csv')
#   
#   raw_files = list.files('data/raw/',full.names = T)
#   
#   all_data = list()
#   
#   for(i in raw_files){
#     # temperature
#     temp = gsub('.csv','',strsplit(i,' ')[[1]][3])
#     if(temp == '36C') next
#     # read in data
#     dat = read.csv(i)
#     # ods
#     ods = get_data(dat,50,51,195)
#     # ec fluorescence (510_540)
#     ecf = get_data(dat,199,200,344)
#     # pp fluorescence (585_615)
#     ppf = get_data(dat,348,349,493)
#     
#     # flatten out stationary phase
#     ec_reps = c('EC1','EC2','EC3')
#     for(j in ec_reps){
#       max_fl = trunc_flatten$fl_flat[trunc_flatten$Temp == temp & trunc_flatten$Replicate == j]
#       if(is.na(max_fl)) next
#       end_exp = min(ecf$Time[ecf[,j] > max_fl])
#       ecf[,j][ecf$Time>=end_exp] = max_fl
#       
#     }
#     pp_reps = c('PP1','PP2','PP3')
#     for(j in pp_reps){
#       max_fl = trunc_flatten$fl_flat[trunc_flatten$Temp == temp & trunc_flatten$Replicate == j]
#       if(is.na(max_fl)) next
#       end_exp = min(ppf$Time[ppf[,j] > max_fl])
#       ppf[,j][ppf$Time>=end_exp] = max_fl
#     }
#     co_reps = c('ECPP1','ECPP2','ECPP3')
#     for(j in co_reps){
#       r = gsub('ECPP','',j)
#       ecco = paste0('ECco',r)
#       ppco = paste0('PPco',r)
#       max_fl_ec = trunc_flatten$fl_flat[trunc_flatten$Temp == temp & trunc_flatten$Replicate == ecco]
#       if(!is.na(max_fl_ec)){
#         end_exp_ec = min(ecf$Time[ecf[,j] > max_fl_ec])
#       }else{
#         max_fl_ec = ecf[nrow(ecf),j]
#         end_exp_ec = ecf$Time[nrow(ecf)]
#       }
#       max_fl_pp = trunc_flatten$fl_flat[trunc_flatten$Temp == temp & trunc_flatten$Replicate == ppco]
#       if(!is.na(max_fl_pp)){
#         end_exp_pp = min(ppf$Time[ppf[,j] > max_fl_pp])
#       }else{
#         max_fl_pp = ppf[nrow(ppf),j]
#         end_exp_pp = ppf$Time[nrow(ppf)]
#       }
#       end_exp = min(end_exp_ec,end_exp_pp)
#       #end_exp = trunc_flatten$end_exp[trunc_flatten$Temp==temp & trunc_flatten$Replicate %in% ecco]
#       ecf[,j][ecf$Time>=end_exp] = ecf[,j][ecf$Time==ecf$Time[ecf$Time >= end_exp][1]]
#       ppf[,j][ppf$Time>=end_exp] = ppf[,j][ppf$Time==ppf$Time[ppf$Time >= end_exp][1]]
#     }
#     # for(j in co_reps){
#     #   r = gsub('ECPP','',j)
#     #   ecco = paste0('ECco',r)
#     #   ppco = paste0('PPco',r)
#     #   max_fl = trunc_flatten$fl_flat[trunc_flatten$Temp == temp & trunc_flatten$Replicate %in% c(ecco,ppco)]
#     #   if(is.na(max_fl)) next
#     #   end_exp = min(ppf$Time[ppf[,j] > max_fl])
#     #   ppf[,j][ppf$Time>=end_exp] = max_fl
#     # }
#     
#     # initialize dataframe
#     all_in_ecf = data.frame(matrix(NA,nrow=nrow(ods),ncol=ncol(ods)+3))
#     colnames(all_in_ecf) = c('Time','Temp','EC1','EC2','EC3','ECco1','ECco2','ECco3','PPco1','PPco2','PPco3','PP1','PP2','PP3')
#     
#     # convert PP fluorescence units to EC fluorescence units (monoculture)
#     
#       max_ecod_mono = mean(unlist(ods[nrow(ods),c('EC1','EC2','EC3')]))
#       max_ppod_mono = mean(unlist(ods[nrow(ods),c('PP1','PP2','PP3')]))
#       max_ecf_mono = mean(unlist(ecf[nrow(ecf),c('EC1','EC2','EC3')]))
#       max_ppf_mono = mean(unlist(ppf[nrow(ppf),c('PP1','PP2','PP3')]))
#       if(temp == '32C'){
#         max_ecod_mono_32 = max_ecod_mono
#         max_ecf_mono_32 = max_ecf_mono
#         max_ppod_mono_32 = max_ppod_mono
#       }
#       norm_ec = max_ecf_mono_32*max_ecod_mono/max_ecod_mono_32/max_ecf_mono
#       norm_pp = max_ecf_mono_32*max_ppod_mono/max_ecod_mono_32/max_ppf_mono
#     # 
#     # 
#       
#     # add time and temp
#     all_in_ecf[,c('Time','Temp')] = ecf[,c('Time','Temp')]
#     
#     # add EC monoculture and coculture data to dataframe as is
#     all_in_ecf[,c('EC1','EC2','EC3','ECco1','ECco2','ECco3')] = ecf[,c('EC1','EC2','EC3','ECPP1','ECPP2','ECPP3')]#*norm_ec
#     
#     all_in_ecf[,c('PPco1','PPco2','PPco3','PP1','PP2','PP3')] = ppf[,c('ECPP1','ECPP2','ECPP3','PP1','PP2','PP3')]#*norm_pp
#     
#     #print(tail(all_in_ecf))
#     
#     all_in_ecf = all_in_ecf %>% gather(Replicate,Growth,EC1:EC3,PP1:PP3,ECco1:ECco3,PPco1:PPco3)
#     
#     ods = ods %>% gather(Replicate,Growth,EC1:EC3,PP1:PP3,ECPP1:ECPP3)
#     print(ggplot(ods, aes(x=Time,y=Growth,col=Replicate)) + geom_point() + theme_bw() + ggtitle(paste(temp,'ODs')))
#     
#     co_only_ec = filter(all_in_ecf,Replicate %in% c('ECco1','ECco2','ECco3'))
#     co_only_ec_ods = filter(ods,Replicate %in% c('ECPP1','ECPP2','ECPP3'))
#     co_only_ec_ods$Growth = co_only_ec_ods$Growth*20000
#     
#     print(ggplot(data.frame(rbind(co_only_ec,co_only_ec_ods)), aes(x=Time,y=Growth,col=Replicate)) + geom_point() + theme_bw() + ggtitle(temp))
#     
#     co_only_pp = filter(all_in_ecf,Replicate %in% c('PPco1','PPco2','PPco3'))
#     co_only_pp_ods = filter(ods,Replicate %in% c('ECPP1','ECPP2','ECPP3'))
#     co_only_pp_ods$Growth = co_only_pp_ods$Growth*5000
#     
#     print(ggplot(data.frame(rbind(co_only_pp,co_only_pp_ods)), aes(x=Time,y=Growth,col=Replicate)) + geom_point() + theme_bw() + ggtitle(temp))
#     
#     mono_only_ec = filter(all_in_ecf,Replicate %in% c('EC1','EC2','EC3'))
#     mono_only_ec_ods = filter(ods,Replicate %in% c('EC1','EC2','EC3'))
#     mono_only_ec_ods$Growth = mono_only_ec_ods$Growth*20000
#     
#     print(ggplot(data.frame(rbind(mono_only_ec,mono_only_ec_ods)), aes(x=Time,y=Growth,col=Replicate)) + geom_point() + theme_bw() + ggtitle(temp))
#     
#     mono_only_pp = filter(all_in_ecf,Replicate %in% c('PP1','PP2','PP3'))
#     mono_only_pp_ods = filter(ods,Replicate %in% c('PP1','PP2','PP3'))
#     mono_only_pp_ods$Growth = mono_only_pp_ods$Growth*5000
#     
#     print(ggplot(data.frame(rbind(mono_only_pp,mono_only_pp_ods)), aes(x=Time,y=Growth,col=Replicate)) + geom_point() + theme_bw() + ggtitle(temp))
#     
#     #break
#     all_data[[temp]] = all_in_ecf
#   }
#   return(all_data)
# }
