#import data
clickstreams = read.csv("C:/Users/JDPur/OneDrive - Durham University/Project/sample_clickstreams.csv", header = FALSE, sep = ",")


cls <- readClickstreams("C:/Users/JDPur/OneDrive - Durham University/Project/sample_clickstreams.csv", header = TRUE, sep = ",") 
mc <- fitMarkovChain(cls) 
mc2 = fitMarkovChain(cls, order = 2)
startPattern <- new("Pattern", sequence = c("H", "c")) 
predict(mc, startPattern) 
plot(mc)

obs_states = c("A", "B", "C" ,"D", "E", "H", "I", "P", "R", "S", "U", "X", "Y")
#For less data
clssmall = cls[1:20]
mcsmall = fitMarkovChain(clssmall)

clssmallstreams = c()
for(i in 1:20){
  clssmallstreams[i] = c(toString(clssmall[[i]]))
}
clssmallstreams
#chi test, tests for independence?

chiSquareTest(clssmall, mcsmall)

#cluster analysis

clusters = clusterClickstreams(cls, order = 1, centers= 2)

#frequencies
nums = frequencies(clssmall)
colSums(nums)

#heatmaps
maps = function(k){
  par(mfrow = c(1, k))
  clust = clusterClickstreams(cls, order =1, centers = k)
  chains = fitMarkovChains(clust)
  for (x in (1:k)){
    hmPlot(chains[[x]], order = 1)
  }

}

hmPlot(fitMarkovChains(clusters)[[1]], order = 1)

#compute hitting times

hit = function(data, character, not_hit = 0){
  hits = c()
  for (i in 1:length(data)){
    stream = unlist(data[i])
    names(stream) = NULL
    if(length(which(stream == character))){
      hits[i] = min(which(stream == character), na.rm = TRUE)
    }else{
      hits[i] = not_hit
    }
    
  }
  mean(hits, na.rm=TRUE)
}

hitting_times = c()
for (i in obs_states){
  hitting_times = c(hitting_times, hit(cls, i, not_hit = 0))
}
barplot(hitting_times, names.arg = obs_states, col = c("red", "orange", "red", "red", "purple", "blue", "blue", "red", "purple", "blue", "blue", "orange", "orange", "orange"), xlab = "Observation States", ylab = "Average Hitting Time", axes=TRUE,)

#Cluster analysis data

cluster_hits = function(data, n, nh = 0){
  for(i in 1:n){
    clusts = clusterClickstreams(data, order = 1, centers = i)
    hit_E = c()
    hit_X = c()
    for(j in 1:i){
      hit_E[j] = hit(clusts[[1]][[j]], "E", not_hit = nh)
      hit_X[j] = hit(clusts[[1]][[j]], "X", not_hit = nh)
    }
    x = data.frame("Clusters" = 1:i, "E time" = hit_E, "X time" = hit_X)
    cat("\n")
    print(x)
  }
}
#consider boxplots of hitting times

boxplot(hit(clusters[[1]][[1]], "D"), hit(clusters[[1]][[2]], "D"))
