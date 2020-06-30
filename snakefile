# snakefile to run everything

import numpy as np

temps = ['32C','34C','36C','38C','40C']
aij = np.round(np.arange(-5,5.1,0.5),2)
aji = np.round(np.arange(-5,5.1,0.5),2)

rule all:
  input:
    expand('results/figures/growth/{temp}_growth.pdf',temp=temps),
    expand('results/monoculture_parameters.csv',temp=temps) ,
    expand('results/coculture/fit_coculture_aij_{aij}_aji_{aji}.csv',aij=aij,aji=aji),
    'results/all_results.csv',
    'results/figures/error/error.pdf'

rule install_github_packages:
  output:
    'results/installed_packages.csv'
  script:
    'lib/install_github_packages.R'

rule process_raw_data:
  input:
    'results/installed_packages.csv',
    expand('data/raw/ECPP_{temp}.csv',temp=temps)
  output:
    expand('data/normalized/{temp}.csv',temp=temps),
    expand('results/figures/growth/{temp}_growth.pdf',temp=temps)
  script:
    'lib/process_raw_data_script.R'

rule fit_monoculture_curves:
  input:
    dat=expand('data/normalized/{temp}.csv',temp=temps)
  output:
    'results/monoculture_parameters.csv'
  script:
    'lib/fit_monoculture_curves.R'

rule fit_coculture_curves:
  input:
    params='results/monoculture_parameters.csv',
    dat=expand('data/normalized/{temp}.csv',temp=temps)
  output:
    'results/coculture/fit_coculture_aij_{aij}_aji_{aji}.csv',
  script:
    'lib/fit_coculture_curves.R'

rule cat_results:
  input:
    expand('results/coculture/fit_coculture_aij_{aij}_aji_{aji}.csv',aij=aij,aji=aji)
  output:
    'results/all_results.csv'
  shell:
    'cat results/coculture/fit_coculture_aij_* | head -n1 > results/all_results.csv; cat results/coculture/* | grep -v temp >> results/all_results.csv'

rule plot_errors:
  input:
    expand('results/coculture/fit_coculture_aij_{aij}_aji_{aji}.csv',aij=aij,aji=aji)
  output:
    'results/figures/error/error.pdf'
  script:
    'lib/plot_errors.R'
