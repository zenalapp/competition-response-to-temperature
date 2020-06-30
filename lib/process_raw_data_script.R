# save normalized data so only have to read it in

# load function
source('lib/process_raw_data.R')

# read in raw data
raw_dat = read_raw_data(dir = 'data/raw/')
# load data about where to truncate (for start of exponential phase) and flatten (start of stationary phase)
trunc_flatten = read.csv('data/trunc_flatten.csv')
# flatten stationary
flattened = flatten_stationary(raw_dat,trunc_flatten)
# normalize data
normalized = normalize_fluorescence(flattened)
# plot raw data
plot_plots(raw_dat,normalized)

# save normalized data
temps = c('32C','34C','36C','38C','40C')
for(i in temps){
  write_csv(x = normalized[[i]], path = paste0('data/normalized/',i,'.csv'))
}
