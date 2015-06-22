genhurst <- function(S,q=1,maxT=19)
{
  L <- length(S)
  H <- c() 
  k = 0
  mcord <- c()   
  lasttmax <- 0
  for(Tmax in 5:maxT)
  {
    k <- k + 1
    x <- 1:Tmax
    for(tt in (lasttmax+1):Tmax)
    {
      dV <- S[seq(tt+1,L,tt)] - S[seq(tt+1,L,tt)-tt]
      VV <- S[seq(tt+1,L+tt,tt)-tt]
      N <- length(dV) + 1
      X <- 1:N
      X <- as.numeric(X)
      Y <- VV
      mx <- sum(X)/N
      #print(paste(N,tt))
      SSxx <- sum(X*X)-N*mx*mx
      my <- sum(Y)/N
      SSxy <- sum(X*Y) - N * mx * my
      cc1 <- SSxy/SSxx
      cc2 <- my - cc1*mx
      ddvd <- dV - cc1
      #print(paste(length(VV),length(dV), length(cc1),N,length(cc2),L))
      vvvd <- VV - cc1*(1:N) - cc2
      mcord <- c(mcord,mean(abs(ddvd)^q)/mean(abs(vvvd)^q))
    }
    lasttmax <- Tmax
    mx <- mean(log10(x))
    SSxx <- sum(log10(x)^2) - Tmax*mx*mx
    my <- mean(log10(mcord))
    SSxy <- sum(log10(x) * log10(mcord)) - Tmax * mx * my
    H[k] <- SSxy/SSxx
  }
  mH <- mean(H)/q
  return(mH)
}

