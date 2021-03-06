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
