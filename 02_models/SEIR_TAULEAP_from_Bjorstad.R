emat2 = matrix(c(1, 0, 0, 0,
                 -1, 0, 0, 0,
                 -1, 1, 0, 0,
                 0, -1, 0, 0,
                 0, -1, 1, 0,
                 0, 0, -1, 0,
                 0, 0, -1, 1,
                 0, 0, 0, -1), ncol = 4, byrow = TRUE)

rlist2 = c(quote(mu * (S + E + I + R)),
           quote(mu * S),
           quote(beta * S * I/(S + E + I + R)),
           quote(mu * E),
           quote(sigma * E),
           quote(mu * I),
           quote(gamma * I),
           quote(mu * R))


tau = function(rateqs, eventmatrix, parameters, initialvals, deltaT, endT){
  time = seq(0, endT, by = deltaT) 
  res = data.frame(matrix(NA, ncol = length(initialvals) + 1, nrow = length(time)))
  res[, 1] = time 
  names(res) = c("time", names(inits)) 
  res[1,]= c(0, inits)
  for (i in 1:(length(time) - 1)) { 
    # calculate overall rates
    rat = sapply(rateqs, eval, as.list(c(parameters, res[i, ]))) 
    evts = rpois(1, sum(rat) * deltaT) 
    if (evts > 0) { 
      # draw events 
      whichevent = sample(1:nrow(eventmatrix), evts, prob = rat, replace = TRUE)
      mt = rbind(eventmatrix[whichevent, ], t(matrix(res[i, -1]))) 
      mt = matrix(as.numeric(mt), ncol = ncol(mt)) 
      # update states
      res[i + 1, -1] = apply(mt, 2, sum)
      res[i + 1, ][res[i + 1,]< 0] = 0 
    } else { 
        # if no events in deltaT
      res[i + 1, -1] = res[i, -1] 
    }
  } 
return(res) 
}


paras = c(mu = 1, beta = 1000, sigma = 365/8, gamma = 365/5) 
inits = c(S = 999, E = 0, I = 1, R = 0) 
sim2 = tau(rlist2, emat2, paras, inits, 1/365, 2)

matplot(sim2[, 1], sim2[, 2:5], type = "l", log = "y", ylab = "Numbers", xlab = "Time") 
legend("bottomright", c("S", "E", "I", "R"), lty = c(1, 1, 1, 1), col = c(1, 2, 3, 4))
    
    
    
    
    
    
    