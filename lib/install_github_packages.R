library(devtools)
library(tidyverse)

devtools::install_github("sprouffske/growthcurver")

write_csv(data.frame(package='growthcurver'),'results/installed_packages.csv')
