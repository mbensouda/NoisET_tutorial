# Inferring responding immune clonotypes
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('clone.R')
# as noisets is not properly called in R so we need to source it one by one.
#source_python('./r-noisets/lib/python3.7/site-packages/noisets/noisettes.py')

#3.1 - Loading data
#For this analysis we will focus on the first time point, day 15, and day 45.
#This tables are imported below using a noiset utility.

#Notice that by default noiset looks for expansion, while here we look for contraction.
#A simple trick for fixing this is just to invert the time points.
# Data information 
path <- file.path(directory, 'NoisET_tutorial/data/')
filename1 <- 'MP_45_F1_short.zip' # first biological replicate
filename2 <- 'MP_15_F2_short.zip' # second biological replicate

#colnames that will change if you work with a different data-set
colnames1 <- c('Clone fraction','Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3') 
colnames2 <- c('Clone fraction','Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3')
# check 
MP_contraction <- ns$Data_Process(path, filename1, filename2, colnames1,  colnames2)
print(paste0("First Filename is : " , MP_contraction$filename1))
print(paste0("Second Filename is : ",  MP_contraction$filename2))
print(paste0("Name of the columns of first file are : ", MP_contraction$colnames1))
print(paste0("Name of the columns of second file are : ", MP_contraction$colnames2))

# Effectively importing the data
ndf <- MP_contraction$import_data()
n <- ndf[[1]]
df_contraction <- ndf[[2]]

print(head(df_contraction, 20))
X = linspace(0,1, 1000)
ggplot() + geom_point(data = df_contraction, aes(x = Clone_fraction_2, y = Clone_fraction_1), color = 'blue', alpha = 0.1) +
  scale_y_continuous(trans='log', labels = scales::label_scientific(1)) + 
  scale_x_continuous(trans='log', labels = scales::label_scientific(1)) + 
  labs(x = 'clone frequency day 15', y = 'clone frequency day 45') +
  coord_cartesian( xlim = c(2.5*10^(-7), 10^(-1)), ylim =  c(2.5*10^(-7), 10^(-1))) +
  geom_line(data = as.data.frame(X), aes(x = X, y = X), linetype = "dashed")


###Task 3.1 Plot in the same figure scatter plots showing the variability between
#replicates samples from the same day, and samples taken at two different time points
#(here day 15 and day 45).

### 3.2 Detection of contracting clones from day 15 to day 45
#Parameters of the noise model learnt in part II.
paras = c(-1.97822857,   1.25456411,   1.04465803, -10.14630235)
# Negative Binomial Sampling Noise Model
noise_model <- 1 

#Below the noiset utility for identifying expanded clonotypes is run. Two parameters
#can be tuned : pval_threshold is a threshold that controls the permissiveness of the
#method (playing somehow the role of the statistical significance) and smed_threshold
#instead can control the weight to give to large colonotypes in being significant.

expansion <- ns$Expansion_Model() # Creating an object for which the associated methods are linked to the contraction/expansion
pval_threshold <- 0.05  #Parameters to play with
smed_threshold <- 0 # Parameters to play with

outpath <- 'contracted_clones' # name of the file Chose what you want 

#Learn the contraction/expansion model + compute different statistics of the log-fold change variable + detect contracting clones
#This part should take approximatively 10 minutes
expansion$expansion_table(outpath, paras, paras, df_contraction, noise_model + 1, pval_threshold, smed_threshold)

#Read the results 
#Board interlude to explain the different statistics 
table_expansion <- read.csv('contracted_clonestop_expanded.csv', sep = '\t')
print(head(table_expansion, 20))

#We can visualise the contracted clonotypes identified by noiset as the red points 
#in the scatter plot below.

ggplot() + geom_point(data = df_contraction, aes(x = Clone_fraction_2, y = Clone_fraction_1), color = 'blue', alpha = 0.1) +
  geom_point(data = table_expansion, aes(x = X.f_2., y = X.f_1.), color = 'red', alpha = 0.2) +
  scale_y_continuous(trans='log', labels = scales::label_scientific(1)) + 
  scale_x_continuous(trans='log', labels = scales::label_scientific(1)) + 
  labs(x = 'clone frequency day 15 after COVID infection', y = 'clone frequency day 45 after COVID infection') +
  coord_cartesian( xlim = c(2.5*10^(-7), 10^(-1)), ylim =  c(2.5*10^(-7), 10^(-1)))

###Task 3.2 Try another value of smed_threshold (name the outpath differently doing so)
#and compare the results in a same figure showing scatter plots for the two different
#thresholds as in the figure above.

###Task 3.3 Do the same analysis for contracted clones between day 15 and day 80 and 
#show the results in a figure. (Compare the number of detected clones for both
#analysis : contracted clones between day 15/day 45, and contracted clones between
#day 15 and day 80)

###3.3 Validating the identified clonotypes
###In this part, we want to check that clonotypes detected by NoisET are actually
#biologically relevant. To do so, we can verify that there is a significant overlap
#between them and the open-access database of TCR from : 
#https://clients.adaptivebiotech.com/pub/covid-2020 known to interact with SARS-Cov2
#antigen. (Here the paper describing this data
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7418734/)

#Let's import the table!

import_mira_clones <- function(path){
  mira_clones <- read.csv(path)['TCR.BioIdentity']
  uniques <- unique(mira_clones['TCR.BioIdentity'])
  mira_clones <- data.frame(list(uniques))
  mira_clones_split <- str_split_fixed(mira_clones$TCR.BioIdentity, "\\+", 3)
  mira_clones['aaCDR3'] <- mira_clones_split[,1]
  mira_clones['Vgene'] <- mira_clones_split[,2]
  mira_clones['Jgene'] <- mira_clones_split[,3]
  return(subset(mira_clones, select = -TCR.BioIdentity))
}
mira_cl <- import_mira_clones('./NoisET_tutorial/data/MIRA_clonotypes.csv')$aaCDR3
print(paste('Number of clonotypes in the MIRA database:', nrow(mira_cl)))

#We are going to use the aminoacid sequences instead of the nucleotide ones
#(it's faster to compare amino-acid than nucleotide strings). For counting the
#overlap between two lists of sequences, we are going to use an external software.
#Below its installation. (This software enables to compare in a reasonable way two
#lists of aminoacid sequences).
#Now we can count how many clonotypes of our noiset list are also present in the mira list.
build_shared_clones <- function(mira_clones, clones, h_dist){
  # Finding matches between amino-acid sequences given Hamming distance
  trie <- atriegc$TrieAA()
  trie$insert_list(mira_clones)        
  shared_clones <- trie$shared_elements(np$asarray(clones), as.integer(h_dist))
  # Building dictionary {candidate_clonotype : list_of_mira_clnotype_matches}
  shared_clones_d <- hash()
  for(cl in 1:length(shared_clones)){
    cl_mira = shared_clones[cl]
    if(!(cl %in% keys(shared_clones_d))){
      shared_clones_d[cl] = c(cl_mira)
    } else{
      shared_clones_d[cl] <- c(shared_clones_d[cl], cl_mira)
    }
  }
  return(shared_clones_d)
}
noiset_cl <- table_expansion$CDR3_AA
shared_clones <- build_shared_clones(mira_cl, noiset_cl, 0)
noiset_overlap <- length(shared_clones)/length(noiset_cl)
print(paste('The overlap between the two lists is:', noiset_overlap))

###To asses that our method works, we also need to compare this number with a "dummy" 
#way of selecting interesting clonotypes, i.e. by sampling at random.
#In the following we generate at random a lot of different samples, compute the
#average overlap and its standard deviation. Given these two number we can
#compute a z-score (https://en.wikipedia.org/wiki/Standard_score) to quantify
#how much our method deviates from the dummy one.

all_aa_seqs <- df_contraction$AACDR3

n_trials <- 100
overlaps <- c()
for(i in 1:n_trials){
  sample <- np$random$choice(all_aa_seqs, nrow(noiset_cl), replace=FALSE)
  shared_clones <- build_shared_clones(mira_cl, sample, 0)
  overlaps <- c(overlaps, length(shared_clones)/length(noiset_cl))
} 
z <- abs(noiset_overlap-np$mean(overlaps))/np$std(overlaps)
print(z)
