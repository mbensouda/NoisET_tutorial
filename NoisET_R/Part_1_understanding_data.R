# After having downloaded the repository, the dataset should be located in 'NoisET_tutorial/data/'.
# The data are organised into tables that we are going to import. The name of the file is composed of two capital letters, 
# which are the acronyms of the patient, a number, which is the number of days after the infection and F1/F2, 
# which label the two different sample replicates at the same day for the same patient.
# Below those tables are imported and stored in a dictionary, where the keys are string specifying the patient, th time and the replicate. 
# The values are dataframes, which are read by unzipping the zip files then using read.table to read the data
# It is useful to know also the command list.dirs (after having imported utils) for iterating among all the file names in a directory.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('clone.R')
import_clones <- function(patient, data_folder){
  # """ 
  # Function that imports all the clonotypes of a given patient and stores
  # them in a dictionary. It returns also the list of ordered time points
  # of such tables.
  # """
  times <- c()
  # hash in R is equivalent to dictionary in Python. However, it needs to get double [[]] whenever we call the field
  clones <- hash()
  # Iteration over all the file in the folder
  for(file_name in list.files(data_folder)){
    # If the name before the underscore corresponds to the chosen patient..
    if(unlist(strsplit(file_name, '_'))[1] == patient){
      # Get the name for the files
      fname <- unlist(strsplit(file_name, '.zip'))[1]
      # Import the table
      frame <- read.table(unzip(file.path(data_folder, file_name), paste0(fname, ".txt")), header = T, sep="\t")
      # Store it in a dictionary where the key contains the patient, the time
      # and the replicate.
      clones[[substring(file_name, first = 1, last = nchar(file_name)-10)]] = frame
      
      # Reading the time from the name and storing it
      times <- c(times, (as.integer(unlist(strsplit(file_name, '_'))[2])))
      print(paste('Clonotypes',substring(file_name, first = 1, last = nchar(file_name)-10),'imported'))
    }
  }
  # Sorting the unique times
  times = sort(unique(times))
  return(list(clones, times))
}

data_folder <- file.path(directory, "NoisET_tutorial/data")
patient <- "MP"
cloneResults <- import_clones(patient, data_folder)
clones <- cloneResults[[1]]
times <- cloneResults[[2]]
print(paste0('Time points: [', toString(times), "]"))

# Display one of these tables and check the different fields. Here we focus on the nucleotide string. 
# We will identify a clonotype with that sequence. For this tutorial we will ignore the amino-acid translation.
head(clones[['MP_15_F1']], 10)

#Task 1.1 How many different clonotypes are present in each table?
#To iterate over a dictionary it's useful to use keys(clones) (to get the list of key values) and nrow (to get the size of the frame).

#[LEAVE EMPTY AS EXERCISE]



# 1.2 - Statistics of clonotype abundance
# It is known that the distribution of counts of clonotypes follow a power law with exponent -2. 
# Here we want to check that our data are meaningful and follow this prediction.
# This plot is hard to reproduce in R
########################################
# bins <- logspace(1, 5.5, 20)
# plt$yscale('log')
# plt$xscale('log')
# plt$xlabel('Clone counts')
# plt$ylabel('pdf')
# totalcl <- data.frame()
# 
# # Plotting all the count distribution of all the tables
# for(id_ in keys(clones)){
#   cl <- clones[[id_]] 
#   h <- plt$hist(cl["Clone.count"], log=TRUE, density=TRUE, bins=bins, histtype='step', lw=2)
#   if (length(h[[1]]) < length(h[[2]])){
#     h[[1]] <- c(h[[1]], 0)
#   }
#   cldf <- data.frame(bins = h[[2]], Clone.count = h[[1]])
#   cldf$group <- id_
#   totalcl <- rbind(totalcl, cldf)
# }
# ggplot(totalcl, aes(x = bins, y = Clone.count, group = group)) + #stat_bin(geom = 'step', stat = 'count', breaks = bins) +
#   geom_line(geom = 'step') + 
#   scale_y_continuous(trans='log') + scale_x_continuous(trans='log') +
#   labs(x = 'Clone counts', y = 'pdf') #+ scale_fill_viridis()
# #coord_cartesian(ylim=c(0.5*10^5, 5*10^6))
# cl <- clones[[keys(clones)[2]]]
# x <- x + ggplot(cl, aes(Clone.count)) + stat_bin(geom="step")+
#   scale_y_continuous(trans='log') + scale_x_continuous(trans='log') +
#   labs(x = 'Clone counts', y = 'pdf')
# # As a check, we plot the power law with exponent -2
# plt$plot(bins, 10*bins**(-2), c='black', ls='--', label= paste('$\\propto n^{-2}$'))
# plt$legend(fontsize=14)
# ###########################################

#plot(bins, 10*bins**(-2), type = 'l', log = c('xy'))
#dev.off()
# 1.3 Common clonotypes between replicates and between time-points
# Here we address the following question: are the clonotypes persistent through time points (i.e. the same sequence appears at multiple times) 
# and through replicates (i.e. how many common sequences the two replicates have)?
# Let's focus on similarity between replicates first. A simple quantification of this is to count common sequences. Do this in the task below.
# Task: how many clonotypes are in common between the two replicates of each time point? Compute also the Jaccard similarity bewteen the replicates, 
# i.e. number of clonotypes in the intersection divided by the number in the union. (https://en.wikipedia.org/wiki/Jaccard_index)
# To count this number it is useful to use python sets and the intersection and union commands (https://www.w3schools.com/python/ref_set_intersection.asp)

#[LEAVE EMPTY AS EXERCISE]


#For quantifying persistence through time points it's useful 
#to merge the replicates in a single table, where the counts of 
#the common clonotypes are summed together. For this purpose two 
#important functions are used: inner_join 

# Creating the dataframes for the merged replicates. After this operation the 
# clones_merged dictionary contains the the merged table of the first and second
# replicate. The indexes are the same as before without the F1/2 label.

clones_merged = hash()

# Iteration over the times
for(t in times){
  # Building the ids correponding at 1st and 2nd replicate at given time point
  id_F1 <- paste0(patient, '_', toString(t), '_F1')
  id_F2 <- paste0(patient, '_', toString(t), '_F2')
  # Below all the rows of one table are appended to the rows of the other
  merged_replicates = full_join(clones[[id_F1]], clones[[id_F2]])
  # But there are common clonotypes that now appear in two different rows 
  # (one for the first and one for the second replicate)! 
  # Below we collapse those common sequences and the counts of the two are summed 
  merged_replicates = merged_replicates %>% group_by(N..Seq..CDR3) %>% summarise(Clone.count = sum(Clone.count)) %>% as.data.frame()
  # The merged table is then added to the dictionary
  clones_merged[[paste0(patient, '_', toString(t))]] = merged_replicates
}

#Now that you have a big merged table for each time point, 
#you can repeat the task above above computing common sequences,
#but for adjacent time points.
#Task1.2 How many common clonotypes there are between adjacent time
#points? And what is the Jaccard index between them?
#Use the merged-replicate dataframes for this analysis.

#[LEAVE EMPTY AS EXERCISE]


#You can notice that the Jaccard distances between time points
#and between replicates are surplisingly similar! It's like the 
#immune system is not changing much in time. The reason behind 
#this is that the change in response to the disease involves only 
#few clonotypes, and the bulk statistics is covering this signal. 
#This is why we'll need to use quite refined tools, like Noiset,
#for identifying the immune response.
#1.4 Persistent clonotypes between time-points
#Let's now focus more on persistence of clonotypes, specifically 
#we want to count, for each unique clonotype, in how many time
#points it is present. Let's call this number its time occurrence.
#There are multiple ways of computing this numbers in python. 
#Here we use the function unique of numpy with the option return
#_counts=True 

#A list of all the clonotypes appearing in all the time points is created
# Note that if one clonotype is present in 2 or more points, it will be repeated
# twice in the list.
all_clones <- data.frame()
for(id_ in keys(clones_merged)){
  cl <- clones_merged[[id_]]
  all_clones <- rbind(all_clones, cl['N..Seq..CDR3'])
}

# The following function returns the list of unique clonotypes and the number of
# repetitions for each of them. 
# Note that the number of repetitions is exactly the time occurrence
unique_clones_table <- as.data.frame(table(all_clones))
unique_clones <- array(unique_clones_table['all_clones'])
time_occurrence <- array(unique_clones_table['Freq'])
time_occurrence_table <- data.frame(time_occurrence = time_occurrence)
ggplot(time_occurrence_table,aes(time_occurrence)) + geom_histogram(stat = 'count', binwidth = 10) +
  scale_y_continuous(trans='log', breaks = c(10^5, 10^6)) +
  labs(x = 'Time Occurrence', y = 'Counts') +
  coord_cartesian(ylim=c(0.5*10^5, 5*10^6))

#1.4 Trajectories of clonotypes
#Here we want to plot the time trajectories of the 50 most abundant
#clonotypes of each time point. More specifically, we want to consider
#the union of the top clonotypes in each of the 5 points and see how
#their count vary in time.

get_top_clones_set <- function(n_top_clones){
  # """
  # This returns the union of the n_top_clones of each time points.
  # """
  ###### QUESTION HERE?
  top_clones <- c()
  for (id_ in keys(clones_merged)){
    cl <- clones_merged[[id_]]
    id_ <- paste0(patient, '_', toString(t)) # where does t come from?
    top_clones_at_time <- levels(factor(cl[order(-cl$Clone.count),'N..Seq..CDR3'][1:n_top_clones]))
    top_clones <- union(top_clones, top_clones_at_time)
  }
  return(unique(top_clones))
}
  

build_traj_frame <- function(top_clones_set){
  # """
  # This builds a dataframe containing the count at all the time points for each 
  # of the clonotypes specified in top_clones_set.
  # The dataframe has also a field that contains the cumulative count.
  # """
  # The trajectory dataframe is initialised with indexes as the clonotypes in
  # top_clones_set
  traj_frame = data.frame(row.names = top_clones_set)
  traj_frame['Clone cumul count'] <- 0
  
  for(id_ in keys(clones_merged)){ 
    cl <- clones_merged[[id_]]
    # Getting the time from the index of clones_merged
    t <- unlist(strsplit(id_, '_'))[2] # t is already a string
    # Selecting the clonotypes that are both in the frame at the given time 
    # point and in the list of top_clones_set
    top_clones_at_time <- intersect(top_clones_set, levels(factor(unlist(cl['N..Seq..CDR3']))))
    # Creating a sub-dataframe containing only the clone in top_clones_at_time
    row.names(cl) <- cl$N..Seq..CDR3
    cl <- cl[-1] # delete the N..Seq..CDR3 column as it already became the index column
    #clones_at_time <- cl[top_clones_at_time, ]
    # Creating a new column in the trajectory frames for the counts at that time
    traj_frame[top_clones_at_time, paste0('t', t)] <- cl[top_clones_at_time, 'Clone.count']
    # The clonotypes not present at that time are NaN. Below we convert NaN in 0s
    traj_frame[is.na(traj_frame)] <- 0
    # The cumulative count for each clonotype is updated
    traj_frame['Clone cumul count'] <- traj_frame['Clone cumul count'] + traj_frame[paste0('t', t)]
  }
  return(traj_frame)
}

top_clones <- get_top_clones_set(50)
traj_frame <- build_traj_frame(top_clones)

#plots
log_counts <- log10(traj_frame['Clone cumul count'])
max_log_count <- max(log_counts)
min_log_count <- min(log_counts)
traj_plot_data <- data.frame() #an initialization for plot
for(i in row.names(traj_frame)){
  row <- traj_frame[i,]
  traj <- as.vector(t(row[c('t15', 't30', 't37', 't45', 't80')])[,1])
  log_count <- log10(unname(row[i,'Clone cumul count']))
  norm_log_count <- (log_count-min_log_count)/(max_log_count-min_log_count)
  traj_tb <- data.frame(times)
  traj_tb$traj <- as.numeric(traj) + 1
  traj_tb$group <- i
  traj_tb$normalized <- norm_log_count
  traj_tb$logCounts <- log_count
  traj_plot_data <- rbind(traj_plot_data, traj_tb)
}
ggplot(traj_plot_data, aes(x = times, y = traj, group = group, color = normalized)) +
  geom_line() + scale_y_continuous(trans='log') + 
  scale_color_viridis(name = 'Log10 cumulative count', guide = "colourbar", breaks = c(0, 1), labels = c(round(min(traj_plot_data$logCounts),1), round(max(traj_plot_data$logCounts),1))) +
  labs(x = 'Time', y = 'Counts')


#What do you notice about the behavior of this largest clonotypes? Are they "interesting"?
#Task1.3 Explore how this behavior changes if you consider less common clonotypes
#We provide a function below for selecting the top abundant clonotypes between a minimum and a maximum position
get_range_top_clones_set <- function(min_n_top_clones, max_n_top_clones){
  # """
  # This returns the union of the clonotypes between a minumum and a maximum 
  # position in the ranking of counts of each time points.
  # """
  top_clones = c()
  for(id_ in keys(clones_merged)){
    cl <- clones_merged[[id_]]
    id_ <- paste0(patient, '_', toString(t))
    top_clones_at_time <- cl.sort_values('Clone count', ascending=False)[min_n_top_clones:max_n_top_clones]
    top_clones_at_time <- levels(factor(cl[order(-cl$Clone.count),'N..Seq..CDR3'][min_n_top_clones:max_n_top_clones]))
    top_clones <- union(top_clones, top_clones_at_time)
  }
    return(top_clones)
    
}
  
#1.5 Making some sense of these trajectories: PCA
#We can try to address the observation that there are some typical different
#behoviors of clonotype trajectories in a more systematic way. Specifically, 
#can we identify clusters of trajectories that follow a similar patterns?
#To this end we want to reduce the dimensions of the trajectories through PCA.
#First, we need to select only the most abundant clonotypes at any time point in the same way as before. In other words, for each time point we want to isolate the top 1000 clonotypes and we consider the union of them. We also want to process the trajectories of these clonotypes by normalizing them by the maximum value of each. This has been done because we are not interested in the overall abundance, but only on the trend of the trajectory.

# Getting the top 1000 clonotypes ate each time point
top_clones <- get_top_clones_set(1000)
# Building trajectory dataframe
traj_frame <- build_traj_frame(top_clones)
# Converting it in a numpy matrix
traj_matrix <- traj_frame[,c('t15', 't30', 't37', 't45', 't80')] %>%
  rowwise() %>%
  mutate(maxval = max(across())) %>%
  as.data.frame()
# Normalize each trajectory by its maximum
norm_traj_matrix <- traj_matrix/traj_matrix$maxval
norm_traj_matrix <- unname(subset(norm_traj_matrix, select = -maxval))

#PCA can be now performed on this set of trajectories. After PCA, we do a 
#hierarchical clustering in the space of the two principal components.
pca <- skldec$PCA(n_components = as.integer(2))
pca <- pca$fit(t(norm_traj_matrix))
clustering <- sklcluster$AgglomerativeClustering(n_clusters= as.integer(3))
clustering <- clustering$fit(t(pca$components_))
clusterdf <- data.frame(x = double(), y = double(), cluster = character())
for (c_ind in 1:clustering$n_clusters){
  x <- pca$components_[1,][clustering$labels_ == (c_ind - 1)]
  y <- pca$components_[2,][clustering$labels_ == (c_ind - 1)]
  tempcluster <- data.frame(list(x, y))
  colnames(tempcluster) <- c('x', 'y')
  tempcluster$cluster <- as.character(c_ind)
  clusterdf <- rbind(clusterdf, tempcluster)
}
clusterdf$cluster <- as.factor(clusterdf$cluster)
ggplot(clusterdf, aes(x = x, y = y, group = cluster, color = cluster)) +
  geom_point(alpha = 0.2) + labs(title = paste0('PCA components (', nrow(norm_traj_matrix), ' trajs)'),
                                  x = paste0('First component (expl var: ', round(pca$explained_variance_ratio_[1],2), ')'),
                                 y = paste0('Second component (expl var: ', round(pca$explained_variance_ratio_[2],2), ')')) +
  theme(legend.position = 'None')

n_cl <- clustering$n_clusters
plotList <- list()
for(cl in 1:n_cl){
  trajs <- norm_traj_matrix[clustering$labels_ == (cl - 1),]
  # separate each traj and components
  traj_combined <- data.frame()
  for(i in 1:nrow(trajs)){
    temptraj <- data.frame(times)
    temptraj$traj <- 0
    temptraj$traj <- as.list(Filter(Negate(is.null), trajs[i,]))
    temptraj$group <- as.character(i)
    traj_combined <- rbind(traj_combined, temptraj)
  }
  traj_combined$viridis <- cl/n_cl
  traj_combined$traj <- as.double(traj_combined$traj)
  p1 <- ggplot(traj_combined, aes(x = times, y = traj, group = group, color = viridis)) + 
    geom_line(alpha = 0.2) + scale_color_viridis(begin = 0, end = cl/n_cl) + 
    labs(x = 'Time', y = 'Normalized counts') + xlim(15,80) + 
    theme(legend.position = 'none')
  
  meanVal <- aggregate(traj_combined[, 'traj'], list(traj_combined$times), mean)
  sdVal <- aggregate(traj_combined[, 'traj'], list(traj_combined$times), sd)
  for(timei in times){
    traj_combined[which(traj_combined$times == timei), 'lower'] <- 
      meanVal[meanVal$Group.1 == timei, 'x'] - sdVal[sdVal$Group.1 == timei, 'x']
    traj_combined[which(traj_combined$times == timei), 'upper'] <- 
      meanVal[meanVal$Group.1 == timei, 'x'] + sdVal[sdVal$Group.1 == timei, 'x']
    traj_combined[which(traj_combined$times == timei), 'mean'] <- 
      meanVal[meanVal$Group.1 == timei, 'x']
  }
  p2 <- ggplot(traj_combined, aes(x = times, y = traj, color = viridis)) + 
    geom_line(aes(x = times, y = mean)) +
        geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                  position=position_dodge(0.05)) +
    labs(x = 'Time', y = 'Normalized counts') + scale_color_viridis(begin = 0, end = cl/n_cl) +
    theme(legend.position = 'none')
  plotList[[cl]] <- list(p1, p2)
}
#axs[1][cl].fill_between(times, np.quantile(trajs, 0.75, axis=0), np.quantile(trajs, 0.25, axis=0), color=colors[cl])
# print out each plot one by one
print(plotList[[1]][1])
print(plotList[[1]][2])
print(plotList[[2]][1])
print(plotList[[2]][2])
print(plotList[[3]][1])
print(plotList[[3]][2])
# we can use grid.arrange to put all in one plot
