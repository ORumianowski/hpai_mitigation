

library(lhs)
Y <- randomLHS(22, 2) 
Y[,1] <- qunif(Y[,1], 1, 10) 
Y[,2] <- qunif(Y[,2], 20, 30)
