library(seqHMM)
library(ggplot2)
library(markovchain)

#import dataset

clicks = read.table("C:/Users/JDPur/OneDrive - Durham University/Project/sample_clickstreams.csv", header = FALSE, sep = ",", col.names = paste0("V",seq_len(42)), fill = TRUE, stringsAsFactors = FALSE)

clicks$V1 = NULL

#as a seqhmm object
clicks_multi = seqdef(clicks, missing = "")

#with clickstreams padded out
fill = function(data){
  for (i in 1:(dim(data)[1])){
    n = dim(data)[2]
    k = which(data[i,] == "E")
    print(i)
    if(n != k){
      data[i, (k+1):n] = rep("E", n-k)}
  }
  return(data)
}

clicks_same = fill(clicks)
clicks_same = seqdef(clicks_same)


#fitting model using seqhmm

make_model = function(data,M, em_tol = 1e-4, em_restarts = 100, logspace = TRUE, printed = 1){
#start = Sys.time()
  N = length(alphabet(data))
  B_init = simulate_emission_probs(M, N)
  A_init = simulate_transition_probs(M)
  pi_init = simulate_initial_probs(M)
  init_model = build_hmm(data, M, A_init, B_init, pi_init)
  return(fit_model(init_model, em_step = TRUE, local_step = TRUE, control_em = list(restart = list(times = em_restarts),print_level = printed, reltol = em_tol, maxeval = 10000), log_space=logspace))
  #end = Sys.time()
  #return(end-start)
}

#Fitting a model using global optimisation methods
global_model = function(data, M){
  N = length(alphabet(data))
  B_init = simulate_emission_probs(M, N)
  A_init = simulate_transition_probs(M)
  pi_init = simulate_initial_probs(M)
  init_model = build_hmm(data, M, A_init, B_init, pi_init)
  return(fit_model(init_model, em_step = FALSE, local_step = TRUE, log_space=TRUE))
}
model1 = make_model(clicks_same, 4)
model2 = make_model(clicks_multi, 4)

#make a heatmap of the transition matrix
transhm = function(matrix, log=FALSE){
  if(log){
    matrix[which(matrix==0)] = 0.01
    matrix = log(matrix)}
  #points = annotater(matrix > 0.95, matrix)
  shaped = reshape2::melt(as.matrix(matrix))
  
  ggplot(data = shaped, aes(x=from, y=to, fill=value) )+scale_fill_gradient(low = "yellow", high = "red")  +     geom_tile() 
  #+ annotate("text", x = points[1,], y = points[2,], label = rep("X", length(points[1,])))
}

#heatmap function for emmision matrix
emhm = function(matrix, log=FALSE){
  #points = annotater(matrix > 0.95, matrix)
  if(log){
    matrix[which(matrix==0)] = 0.01
    matrix = log(matrix)}
  shaped = reshape2::melt(as.matrix(matrix))
  ggplot(data = shaped, aes(x=state_names, y=symbol_names, fill=(value)) ) +scale_fill_gradient(low = "lightsteelblue1", high = "purple") +     geom_tile() 
  #+ annotate("text", x = points[1,], y = points[2,], label = rep("X", length(points[1,])))
}

#Plot Heatmaps side by side
hms = function(HMM, log_space = FALSE){
require(gridExtra)
plot1 = transhm(HMM$model$transition_probs, log = log_space)
plot2 = emhm(HMM$mode$emission_probs, log = log_space)
plot = grid.arrange(plot1, plot2, ncol=2)
return(plot)}

#create models for different values of M

model2=make_model(clicks_multi,2)
model3=make_model(clicks_multi,3)
model4=make_model(clicks_multi,4)
model5=make_model(clicks_multi,5)
model6=make_model(clicks_multi,6)
model7=make_model(clicks_multi,7)
model8=make_model(clicks_multi,8)
model9=make_model(clicks_multi,9)
models = list(model2,model3, model4, model5, model6, model7, model8, model9)

#create plot of the frequencies of different likelihood for the data 
liks_plot = ggplot(em_likelihoods, aes(x=states, y=model_restarts))+geom_point(shape = 1)+geom_point(data = (final_likelihoods),aes(x = states, y=model_likelihoods), color="red", shape = 25)
plot(liks_plot+geom_point(data = freq_data[3:102,], aes(x = states, y = likelihoods, size = frequencies, fill = frequencies), color = "blue"))


#Test the EM method for different size data
variances = c()
eigens = c()

for (i in (1/(8*rep(0.2, 20))- 0.5)){
  print(i)
  trans_mat = matrix(rbeta(49, i,i), nrow = 7, ncol = 7)
  for(k in 1:2){
  trans_mat[k,] = trans_mat[k,]/sum(trans_mat[k,])}
  variances = c(variances, var(as.vector(trans_mat)))
  eigens = c(eigens, eigen(trans_mat)$values[-1])
}
  
plot(1:20, variances, type = "l")

#function to generate matrices with different uncertainties
structured_matrix = function(pert, N, M){
    M = diag(pert, nrow =M, ncol=N)
    M[M == 0] = (1-pert)/(N-1)
    return(M)
  }
 
#simulate from a markov chain

generate_init_model = function(T, var_A, var_B, var_pi, N, M){
  A = structured_matrix(var_A,M, M)
  B = structured_matrix(var_B,N, M)
  pi = as.vector(structured_matrix(var_pi,M, 1))
  return(build_hmm(simulate_hmm(T, pi, A, B, 10)$observations, M, A, B, pi))
}


#implement HMM distance measure (HSD)
stationary_dist = function(lambda0){
  return(steadyStates(new("markovchain", states = lambda0$state_names, transitionMatrix = lambda0$transition_probs)))
}

HMM_distance = function(lambda0, lambda1){
  N = lambda0$n_symbols
  M = lambda0$n_states
  
  mc0 = new("markovchain", states = lambda0$state_names, transitionMatrix = lambda0$transition_probs)
  mc1 = new("markovchain", states = lambda1$state_names, transitionMatrix = lambda1$transition_probs)
  
  stat0 = steadyStates(mc0)
  stat1 = steadyStates(mc1)
  
  dist0 = c(sum(stat0*lambda0$emission_probs[,1]))
  dist1 = c(sum(stat1*lambda1$emission_probs[,1]))

  for(i in 2:N){
    dist0 = c(dist0, (dist0[i-1] + sum(stat0*lambda0$emission_probs[,i])))
    dist1 = c(dist1, (dist1[i-1] + sum(stat1*lambda1$emission_probs[,i])))
    
  }
  #print(dist0)
  #print(dist1)
  return(sum(abs(dist0 - dist1)))
}

#test the distance measure
m1 = generate_init_model(100, 0.8, 0.8, 0.8, 5, 2)

A = matrix(c(0.79, 0.21, 0.21, 0.79), nrow = 2)
B = matrix(c(0.78, 0.06, 0.06, 0.04, 0.06, 0.04, 0.78, 0.06, 0.06, 0.06), nrow = 2, byrow = TRUE)
pi = c(0.79, 0.21)
m2 = build_hmm(simulate_hmm(100, pi, A, B, 10)$observations, 2, A, B, pi)
#testing our seqhmm method

test_model = generate_init_model(500, 0.9, 0.9, 0.9, 5, 2)
start_time = Sys.time()
start_time = Sys.time()
sim_model = make_model(test_model$observations, 2, em_tol = 1e-4, em_restarts = 500, printed=0)
end_time = Sys.time()
print(end_time-start_tine)


#Generate simulated data
lik_method_df1 = data.frame(Model = numeric(), Observations = numeric(), Restarts = numeric(), Tolerance = numeric(), logLik = numeric(), LikDif = numeric(), Distance = numeric(), Time = numeric())

repeats = c(1, 10, 50, 100, 500)
pertubations = c(0.6, 0.7, 0.8, 0.9)
tolerance = c(10e-2, 10e-3, 10e-4, 10e-5)
Ts = c(1:9, seq(10, 100, by = 10), seq(200, 1000, by = 100))

i = 1
pb = txtProgressBar(min = 0, max = 112, style =3)
init_liks = c()
for(p in pertubations){
  for(l in Ts){
    init_model = generate_init_model(l, p, p, p, 5, 2)
    init_lik = logLik(init_model)
    init_liks = c(init_liks, init_lik)
    #for(r in repeats){
      #for(t in tolerance){
       
        start_time = Sys.time()
        sim_model = global_model(init_model$observations, 2)
        end_time = Sys.time()
        if(sim_model$model$emission_probs[1,1]<sim_model$model$emission_probs[2,1]){
          init_model$emission_probs = init_model$emission_probs[c(2,1),]
        }
        sim_lik = sim_model$logLik
        lik_difference = sim_lik - init_lik
        distance = HMM_distance(init_model, sim_model$model)
        lik_method_df1[nrow(lik_method_df1)+1, ]= c(which(pertubations==p), l, 1, 1e-8, sim_lik, lik_difference, distance, as.double(difftime(end_time, start_time, units = "secs")))
        i= i+1
        setTxtProgressBar(pb,i)
        #}
        
        
      #}
    }
}


#test effect on tollerance
tol_diff_vars = c()

for(i in seq(1, 236, by = 4)){
  tol_diff_vars = c(tol_diff_vars, var(lik_method_df[i:(i+4), 7]))
}
#Checking effect on accuracy

factor_lik_method1 = lik_method_df1
for(i in 1:4){factor_lik_method1[,i] = as.factor(factor_lik_method1[,i])}

#plot of distributions of tolerance and restarts vs accuracy

p1 = ggplot(factor_lik_method1, aes(x=Model, y=Distance, fill = Tolerance))+geom_boxplot()+scale_fill_brewer(palette = "PuRd")
p2 = ggplot(factor_lik_method1, aes(x=Model, y=Distance, fill = Restarts))+geom_boxplot()+scale_fill_brewer(palette = "PuRd")
p3 = ggplot(factor_lik_method1, aes(x=Model, y=LikDif, fill = Tolerance))+geom_boxplot()+scale_fill_brewer(palette = "Blues")
p4 = ggplot(factor_lik_method1, aes(x=Model, y=LikDif, fill = Restarts))+geom_boxplot()+scale_fill_brewer(palette = "Blues")
require(gridExtra)#
grid.arrange(p1, p2,p3, p4, nrow=2)

#plot to assess effect of more observations

ggplot(data = factor_lik_method[which((factor_lik_method$Restarts == 500) & (factor_lik_method$Tolerance == 1e-04 )),], aes(x=Observations, y=Distance, group = Model))+geom_line(aes(color=Model))+geom_point(aes(color=Model))+scale_color_brewer(palette = "Accent")

#Assessing the speed 

#Looking at all possible speeds
ggplot(data = factor_lik_method[which(factor_lik_method$Model==1),], aes(x = Observations, y = Time, group = interaction(Restarts, Tolerance)))+geom_point(aes(color = Restarts, shape = Tolerance))+geom_line(aes(color=Restarts, linetype = Restarts ))

#considering only min and max for each observations
mins = c()
maxs = c()

for(i in c(1:4)){
  for(j in Ts){
    mins = c(mins, ((i-1)*560 + (which(Ts==j) -1)*20) + which.min(factor_lik_method1[which((factor_lik_method1$Model==i)&(factor_lik_method1$Observations==j)),]$Time))
    maxs = c(maxs, ((i-1)*60 + (which(Ts==j) -1)*20)+which.max(factor_lik_method1[which((factor_lik_method1$Model==i)&(factor_lik_method1$Observations==j)),]$Time))
  }
}


#Plotting min and max times
ggplot(data = factor_lik_method[ mins,], aes(x = Observations, y = Time, group = Model))+geom_line(aes(linetype =Model), color='purple')+geom_point(aes(color = Restarts, shape = Tolerance))+scale_color_brewer(palette = "Set1")


       
#Assessing effect on more observations

acc_plots = list()

#for EM + LO
for (i in 1:4){
acc_plots[[i]] = ggplot(data = factor_lik_method1[which((factor_lik_method1$Restarts == 500) & (factor_lik_method1$Tolerance == 1e-04 )&factor_lik_method1$Model==i),], aes(x=Observations, y=Distance,))+geom_line(color="purple")+geom_point(color="magenta")+ggtitle(paste("Model ",i))+scale_x_continuous(breaks = seq(0, 1000, by = 100))+theme_grey()+theme(plot.title = element_text(hjust=0.5))
}
ggarrange(plotlist = acc_plots, nrow=2, ncol=2)

#for LO
for (i in 1:4){
  acc_plots[[i]] = ggplot(data = factor_lik_method1[which((factor_lik_method1$Restarts %in% c(0, 500)) & (factor_lik_method1$Model==i)&(factor_lik_method1$Observations%in%seq(0, 1000, by=10))),], aes(x=Observations, y=Distance, group = Restarts))+geom_line(aes(color=Restarts))+ggtitle(paste("Model ",i))+scale_x_continuous(breaks = seq(0, 1000, by = 100))+theme_grey()+theme(plot.title = element_text(hjust=0.5))
}
ggarrange(plotlist = acc_plots, nrow=2, ncol=2, common.legend = TRUE, legend = "right")

#Assesing number of observations vs speed
full_speed_plots = list()

ggarrange(plotlist = full_speed_plots, ncol=2, nrow=2, common.legend = TRUE, legend = "right")
for(i in 1:4){
full_speed_plots[[i]] = ggplot(data = factor_lik_method1[which(factor_lik_method1$Model==i&factor_lik_method1$Observations%in%c(  100*(1:10))&factor_lik_method1$Restarts%in%c(1,10)),], aes(x = Observations, y = Time, group = interaction(Restarts, Tolerance)))+geom_point(aes(color = Restarts, shape = Tolerance))+geom_line(aes(color=Restarts, linetype = Restarts ))+scale_x_continuous(breaks=c(100*(1:10)))+theme_grey()+theme(plot.title = element_text(hjust=0.5))+ggtitle(paste("Model ",i))
}

ggarrange(plotlist = full_speed_plots, ncol=2, nrow=2, common.legend = TRUE, legend = "right")

#with minimums
min_speed_p1 = ggplot(data = factor_lik_method1[ mins,], aes(x = Observations, y = Time, group = Model))+geom_line(aes(linetype =Model), color='purple')+geom_point(aes(color = Restarts, shape = Tolerance))
#+scale_color_brewer(palette = "Set1")

#minimums per model
min_speed_graphs = list()

for(i in 1:4){
min_speed_graphs[[i]]= ggplot(data = rbind(factor_lik_method1[ mins[(i-1)*28+c(1, 10:28)],]), aes(x = Observations, y = Time))+geom_line(color='purple') +geom_line(data=time_data[28*(i-1)+1:28,], aes(x =Observations, y=Time), color = 'darkgreen', linetype = 'dashed')+geom_point(data=time_data[28*(i-1)+1:28,], aes(x =Observations, y=Time), color = 'darkgreen')+geom_point(aes(color = Restarts, shape = Tolerance))+ggtitle(paste("Model ",i))+scale_x_continuous(breaks=seq(0,1000, by=100))+scale_fill_manual(name="", values=c("darkgreen", "purple"))
}
ggarrange(plotlist = min_speed_graphs, nrow = 2, ncol=2, common.legend=TRUE, legend = "right")

min_speed_graphs2 = ggplot(data = factor_lik_method1[ mins[(2-1)*28+c(1, 10:28)],], aes(x = Observations, y = Time))+geom_line( color='purple')+geom_line(data=time_data[1:28 +28,], aes(x =Observations, y=Time), color = 'darkgreen', linetype = 'dashed')+geom_point(data=time_data[1:28,], aes(x =Observations, y=Time), color = 'darkgreen')+geom_point(aes(color = Restarts, shape = Tolerance))+ggtitle(paste("Model ",2))+scale_x_continuous(breaks=seq(0,1000, by=100))
min_speed_graphs3 = ggplot(data = factor_lik_method1[ mins[(3-1)*28+c(1, 10:28)],], aes(x = Observations, y = Time))+geom_line( color='purple')+geom_line(data=time_data[1:28 + 56,], aes(x =Observations, y=Time), color = 'darkgreen', linetype = 'dashed')+geom_point(data=time_data[1:28+84,], aes(x =Observations, y=Time), color = 'darkgreen')+geom_point(aes(color = Restarts, shape = Tolerance))+ggtitle(paste("Model ",3))+scale_x_continuous(breaks=seq(0,1000, by=100))
min_speed_graphs4 = ggplot(data = factor_lik_method1[ mins[(4-1)*28+c(1, 10:28)],], aes(x = Observations, y = Time))+geom_line( color='purple')+geom_line(data=time_data[1:28,], aes(x =Observations, y=Time), color = 'darkgreen', linetype = 'dashed')+geom_point(data=time_data[1:28,], aes(x =Observations, y=Time), color = 'darkgreen')+geom_point(aes(color = Restarts, shape = Tolerance))+ggtitle(paste("Model ",4))+scale_x_continuous(breaks=seq(0,1000, by=100))

require(gridExtra)
ggarrange( min_speed_graphs1,min_speed_graphs2, min_speed_graphs3, min_speed_graphs4, nrow = 2,ncol = 2, common.legend = TRUE, legend = "right")


#the sequential pruning strategy

seq_pruning = function(observations, kmin, kmax){
  criteria = c()
  k = kmax
  lambda_k = make_model(observations, k, printed=0, em_restarts = 50)
  criteria[k] = BIC(lambda_k$model)
  p_inf = which.min(stationary_dist(lambda_k$model))
  print(stationary_dist(lambda_k$model))
  print(paste("State", p_inf, "with probability",stationary_dist(lambda_k$model)[p_inf], "pruned"))
  #print(lambda_k$model)
  k = k-1
  #print(stationary_dist(lambda_k$model))
  while(k>=kmin){
  new_A = lambda_k$model$transition_probs[-p_inf, -p_inf] + 1e-8
  new_B = lambda_k$model$emission_probs[-p_inf,]  +1e-8
  new_pi = lambda_k$model$initial_probs[-p_inf]+1e-8

  new_A = new_A/rowSums(new_A) 
  new_B = new_B/rowSums(new_B)
  new_pi = new_pi/sum(new_pi) 

  new_lambda = build_hmm(observations, k, new_A, new_B, new_pi)
  
  lambda_k = fit_model(new_lambda, em_step = TRUE, local_step = TRUE, control_em = list(restart = list(times = 50),print_level = 0, reltol = 1e-4, maxeval = 10000), log_space=FALSE)
  criteria[k] = BIC(lambda_k$model)
  p_inf = which.min(stationary_dist(lambda_k$model))
  
  k=k-1
  }
  return(c(which.min(criteria), criteria[which.min(criteria)]))
}

#sequential pruning of our data

M_selection = seq_pruning(clicks_multi, 2, 9)
click_pruning_criteria = c(NA,4819.200, 4827.583, 4766.655, 4830.159, 4874.612, 4935.204,5057.235, 5196.607)

selection_criteria = data.frame(Criteria = character(), M = numeric(), likelihood = numeric(), stringsAsFactors = FALSE)
#plot model selection for different models

for(i in 2:9){
  selection_criteria[nrow(selection_criteria)+1,]= c("HQC", i, HQC(models[[i-1]], 1363))
}

selection_criteria$Criteria = as.factor(selection_criteria$Criteria)
selection_criteria$Value=as.numeric(selection_criteria$Value)

ggplot(data = selection_criteria, aes(x = M, y = Value, group = Criteria))+geom_line(aes(color= Criteria, linetype = Criteria))+geom_point(aes(color = Criteria))+scale_y_continuous()+labs(y = "Criteria Value")+scale_linetype_manual(values = c("solid","solid","solid","dotted","twodash"))
#+scale_color_brewer(palette = "Set1")

total_length = 0
for(i in 1:200){
  total_length = total_length + which(clicks[i,] == "E")
}

#plot the log likelihood of different models

no_em_model = build_hmm(clicks_multi, 13, simulate_transition_probs(13), diag(13), simulate_initial_probs(13))
no_em_model = fit_model(no_em_model, em_step = TRUE, local_step = TRUE, control_em = list(restart = list(times = 50),print_level = 1, reltol = 1e-4, maxeval = 10000), log_space=TRUE)

ggplot(data.frame(Likelihoods = model_likelihoods, M = as.factor(2:9)), aes(x=M, y = Likelihoods, group = 1))+geom_line(color = "purple")+geom_point(color = "red")+geom_hline(yintercept =   model13$logLik, color = "blue", linetype = "dashed")
#+geom_hline(yintercept = no_em_model$logLik, color="red")

#implement the AIC methods
AICc = function(M, T){
  #L = M$logLik
  p = attr(logLik(M$model), "df")
  return(AIC(M$model) + (2*p*(p+1))/(T-p-1))
}

#implement the HQC criteria
HQC = function(M, T){
  L = M$logLik
  p = attr(logLik(M$model), "df")
  return(-2*L + p*log(log(T)))
}




