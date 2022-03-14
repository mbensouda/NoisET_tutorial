# Require R >= 3.6
if (!require('git2r')) install.packages('git2r'); library(git2r)
if (!require('rstudioapi')) install.packages('rstudioapi')
if (!require('hash')) install.packages('hash'); library(hash)
if (!require('utils')) install.packages('utils'); library(utils)
if (!require('pracma')) install.packages('pracma'); library(pracma)
if (!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if (!require('dplyr')) install.packages('dplyr'); library(dplyr)
if (!require('viridis')) install.packages('viridis'); library(viridis)
if (!require('BBmisc')) install.packages('BBmisc'); library(BBmisc)
if(!require('devtools')) install.packages("devtools")
if (!require('stats')) install.packages('stats'); library(stats)
if (!require('gridExtra')) install.packages('gridExtra'); library(gridExtra)
if (!require('data.table')) install.packages('data.table'); library(data.table)
if (!require('reticulate')) install.packages('reticulate'); library(reticulate)
if (!require('scales')) install.packages('scales'); library(scales)
if (!require('stringr')) install.packages('stringr'); library(stringr)
if (!require('grid')) install.packages('grid'); library(grid)
# make sure that NoisET folder is empty
directory <- getwd()
# if you already get data from Github, then comment the line below.
system('git clone https://github.com/mbensouda/NoisET_tutorial.git')
# set virtual environment to run Python packages
tryCatch({
  use_virtualenv("./r-noisets")
}, error=function(cond) {
  print('This is your first time!')
  virtualenv_create("./r-noisets")
  use_virtualenv("./r-noisets")
})
#set python path
Sys.setenv(RETICULATE_PYTHON = "./r-noisets/bin/python")
# need to reload libary reticulate
library(reticulate)
# check if noisets package is already installed and install other required packages
tryCatch({
  ns <- import('noisets.noisettes')
  atriegc <- import('atriegc')
  np <- import('numpy')
  skldec <- import('sklearn.decomposition')
  sklcluster <- import('sklearn.cluster')
  plt <- import('matplotlib.pyplot')
}, error=function(cond) {
  virtualenv_install(envname = "./r-noisets", packages = 'noisets')
  virtualenv_install(envname = "./r-noisets", packages = 'matplotlib')
  virtualenv_install(envname = "./r-noisets", packages = 'pandas')
  virtualenv_install(envname = "./r-noisets", packages = 'numpy')
  virtualenv_install(envname = "./r-noisets", packages = 'seaborn')
  virtualenv_install(envname = "./r-noisets", packages = 'scipy')
  virtualenv_install(envname = "./r-noisets", packages = 'sklearn')
  virtualenv_install(envname = "./r-noisets", packages = './NoisET_tutorial/ATrieGC/')
  ns <- import('noisets.noisettes')
  atriegc <- import('atriegc')
  np <- import('numpy')
  skldec <- import('sklearn.decomposition')
  sklcluster <- import('sklearn.cluster')
  plt <- import('matplotlib.pyplot')
})
