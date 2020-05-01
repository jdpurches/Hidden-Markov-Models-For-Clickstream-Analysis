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
  
  print(logLik(new_lambda))
  lambda_k = fit_model(new_lambda, em_step = TRUE, local_step = TRUE, control_em = list(restart = list(times = 50),print_level = 0, reltol = 1e-4, maxeval = 10000), log_space=FALSE)
  criteria[k] = BIC(lambda_k$model)
  p_inf = which.min(stationary_dist(lambda_k$model))
  print(stationary_dist(lambda_k$model))
  print(paste(lambda_k$model$state_names[p_inf], "with probability",stationary_dist(lambda_k$model)[p_inf], "pruned"))
  
  k=k-1
  }
  print(criteria)
  return(c(which.min(criteria), criteria[which.min(criteria)]))
}
