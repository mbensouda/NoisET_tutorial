# Inferring responding immune clonotypes
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('clone.R')
# as noisets is not properly called in R so we need to source it one by one.
#source_python('./r-noisets/lib/python3.7/site-packages/noisets/noisettes.py')
### 2.1 Loading data
#Below we load the data for two replicates at a given time point (the first).
#We will use a noiset utilities called Data_Process.
path <- file.path(directory, 'NoisET_tutorial/data/')
filename1 <- 'MP_15_F1_short.zip' # first biological replicate
filename2 <- 'MP_15_F2_short.zip' # second biological replicate

#colnames that will change if you work with a different data-set
colnames1 <- c('Clone fraction','Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3') 
colnames2 <- c('Clone fraction','Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3')
# check 
MP_15 <- ns$Data_Process(path, filename1, filename2, colnames1,  colnames2)
print(paste0("First Filename is : " , MP_15$filename1))
print(paste0("Second Filename is : ",  MP_15$filename2))
print(paste0("Name of the columns of first file are : ", MP_15$colnames1))
print(paste0("Name of the columns of second file are : ", MP_15$colnames2))

# Effectively importing the data
ndf <- MP_15$import_data()
n <- ndf[[1]]
df_15 <- ndf[[2]]

# 2.2 Visualizing the replicate noise
# The experimental sampling noise can be visualized in a scatter plot, 
# where each point is a clonotype and the two axis are the frequencies -- 
# normalized abundances -- in the first and second replicate. 
# In a "perfect" experiment -- meaning sequencing the entire TCR repertoire of one 
# individual, one would expect the counting of the same clonotype to be the same. 
# Therefore, deviations from the diagonal shows this uncertainty.
X = linspace(0,1, 1000)
ggplot() + geom_point(data = df_15, aes(x = Clone_fraction_1, y = Clone_fraction_2), color = 'blue', alpha = 0.1) +
  scale_y_continuous(trans='log', labels = scales::label_scientific(1)) + 
  scale_x_continuous(trans='log', labels = scales::label_scientific(1)) + 
  labs(x = 'clone frequency first replicate', y = 'clone frequency second replicate') +
  coord_cartesian( xlim = c(2.5*10^(-7), 10^(-1)), ylim =  c(2.5*10^(-7), 10^(-1))) +
  geom_line(data = as.data.frame(X), aes(x = X, y = X), linetype = "dashed")

#2.3 - Learning sampling noise
#Noiset quantifies sampling noise described previousely by learning a model.

#There are 3 possible type of models to choose, which at increasing level of 
#complexity are: the Poisson model (2 parameters), the Negative Binomial model 
#(3 parameters), and the Negative Binomial + Poisson model (5 parameters).

#Details about the distributions and the meaning of each parameters are detailed
#in methods section of READme document of the github repository. On this part of
#the tutorial, we will focus on the negative binomial model, which provides reliable
#results in a resonable computational time.

#The Noiset function starts from a given set of initial parameters, and returns 
#the same number of parameters as a result of the optimizaition process.

noise_model <- 1 # Negative Binomial 

#other models :
# 0 : NB + Poisson
# 1 : NB
# 2 : Poisson

# Suggested initial parameters (they are very close to the optimal one 
# to speed up the computation)
init_paras_arr = list( 
  c(-2.0, 1.54, 1.23, 6.65, -9.71),
  c(-2.0, 1.26, 1.05, -10.1),
  c(-2.15, -9.47)
)
# python starts from 0 while R starts from 1 so we need to add 1 after noise_model
# 
init_paras <- init_paras_arr[[noise_model + 1]]
null_model <- ns$Noise_Model() 
# the next line may take a while to run
result <- null_model$learn_null_model(df_15, noise_model, init_paras)
# Optimal parameters
result[[1]]['x']
head(result)
#2.4 Generate synthetic data with NoisET
#Noiset allows also to generate synthetic data of clonotype counts 
#which obey to a given noise model. Here we generate a new synthetic dataset
#with the same experimental features (number of reads, number of clonotypes)
#of the original one.

Synthetic <- ns$Generator()
noise_model <- 1 # Negative Binomial Noise model
NreadsI <- sum(df_15['Clone_count_1'])[1] # Total number of reads in the first sample
NreadsII <- sum(df_15['Clone_count_2'])[1] # Total number of reads in the second sample
Nsamp <- nrow(df_15) # total number of clones found in both samples

samples <- Synthetic$gen_synthetic_data_Null(result[[1]]['x'], noise_model + 1, NreadsI, NreadsII, Nsamp)
f_samples <- samples[[1]]
pair_samples <- samples[[2]]

print(head(pair_samples, 10))

#**Task 2.1** 
#Compare the noise through the scatter plot done before between the real data
#and the synthetic data.


