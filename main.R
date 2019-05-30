
# STAT 950 Statistical Machine Learning (Spring 2019)
      
      # The gibbsBLasso function is adapted from 
      # https://cs.gmu.edu/~pwang7/gibbsBLasso.html
      # Any bugs in the Gibbs sampler?
      # One of the bugs is fine with the p<n case
      # but will lead to problems for p>n. 


      rm(list=ls(all=TRUE))
      # library(pscl) # rigamma
      install.packages(c("mnormt","VGAM","miscTools"))
      library(mnormt) # rmnorm
      library(VGAM) # rinv.gaussian
      library(miscTools) # colMeans
      library(MASS)


gibbsBLasso = function(x, y, max.steps) {
	n <- nrow(x)
	m <- ncol(x)

	XtX <- t(x) %*% x	#Time saving
	xy <- t(x) %*% y

	r <- 1
	delta <- 1.78

	betaSamples <- matrix(0, max.steps, m)
	sigma2Samples <- matrix(1, max.steps,1) 
	invTau2Samples <- matrix(0, max.steps, m)
	lambdaSamples <- matrix(0, max.steps,1)

	beta <- drop(backsolve(XtX + diag(nrow=m), xy))
	residue <- drop(y - x %*% beta)
	sigma2 <- 1 # drop((t(residue) %*% residue) / n) # true sigma2=1
	invTau2 <- 1 / (beta * beta)
	lambda <- m * sqrt(sigma2) / sum(abs(beta))

	k <- 0
	while (k < max.steps) {
		k <- k + 1

		if (k %% 1000 == 0) {
			cat('Iteration:', k, "\r")
		}

		# sample beta
		invD <- diag(invTau2)
		invA <- solve(XtX + invD)
		mean <- invA %*% xy
		varcov <- sigma2 * invA
		beta <- drop(rmnorm(1, mean, varcov))
		betaSamples[k,] <- beta

		# sample sigma2
		shape <- (n+m-1)/2
		residue <- drop(y - x %*% beta)
		scale <- (t(residue) %*% residue + t(beta) %*% invD %*% beta)/2
		sigma2 <- 1/rgamma(1, shape, 1/scale)
		sigma2Samples[k] <- sigma2

		# sample tau2
		muPrime <- sqrt(lambda^2 * sigma2 / beta^2)
		lambdaPrime <- lambda^2
		invTau2 <- rep(0, m)
		for (i in seq(m)) {
			invTau2[i] <- rinv.gaussian(1, muPrime[i], lambdaPrime)
		}
		invTau2Samples[k, ] <- invTau2

		# update lambda
		shape = r + m/2
		scale = delta + sum(1/invTau2)/2
		lambda <- rgamma(1, shape, 1/scale)
		# if (k %% 10 == 0) {
			# low <- k - 9
			# high <- k
			# lambda <- sqrt( 2*m / sum(colMeans(invTau2Samples[low:high, ])) )
		# }
		lambdaSamples[k] <- lambda
	}

      
      # t1 = colMedians(betaSamples[seq(max.steps/2, max.steps, 5), ])
      # t2 = colMedians(sigma2Samples[seq(max.steps/2, max.steps, 5), ])
      t1 = betaSamples
      t2 = sigma2Samples
      t3 = lambdaSamples
	dat = list(t1=t1,t2=t2,t3=t3)
      return(dat)
}


      # real data
      #library(lars)
      #data(diabetes)
      #dim(x)
      #x <- scale(diabetes$x)
      #y <- scale(diabetes$y)



      # simulated data
      n=200; p=10;
      sig = matrix(0,p,p)

      for (i in 1:p)
         {
            for(j in 1:p)
             {
               sig[i,j] = 0.5^abs(i-j)
             }
         }

      x = mvrnorm(n,rep(0,p),sig)
      x = scale(x)
      error = rnorm(n,0,1)
      y = 2.5*x[,1]+3*x[,4]+error

      betas <- gibbsBLasso(x, y, max.steps = 100000)

      max.steps = 100000
      
      # A quick check on the results
      
      apply(betas$t1, 2, median)

      round(apply(betas$t1, 2, median), digits = 2)

      quants = apply(betas$t1, 2, quantile, prob=c(0.025, 0.975))

      abs(sign(quants[1,])+sign(quants[2,])) == 2
      

      # The dth predictor
      # Change d from 1 to 10 to observe the posterior distribution of beta
      # d = 1 and 4 correspond to important predictors
      # Can you obtain sparsity directly based on the posterior estimate?
      
      d = 5 
      
      t1=as.matrix(betas$t1[,d])
      t1 = t1[seq(max.steps/2, max.steps,20),]
      plot(t1)
      median(t1)
      sd(t1)
      quantile(t1,c(0.025,0.5,0.975))
      hist(t1,freq=F,main="The dth predictor")

      my_mode <- function(x){ ## where x is the vector of save iterations of a parameter in stan
           dx <- density(x)
           return(mode=dx$x[which.max(dx$y)])
         }
      my_mode(t1)            ### the posterior mode

      # sigma2
      t2=betas$t2  # true value = 1
      t2 = t2[seq(max.steps/2, max.steps,10),]
      plot(t2)
      median(t2)                      # posterior median estimate
      sd(t2)
      quantile(t2,c(0.025,0.5,0.975)) # 95% posterior credible intervals 
      hist(t2,freq=F,main="The sigma2")
      my_mode(t2)                     ### the posterior mode


     # lambda
     t3=betas$t3  
     t3 = t3[seq(max.steps/2, max.steps,10),]
     plot(t3)
     median(t3)


     #t2=betas$t2
     #t2 = colMedians(t2[seq(max.steps/2, max.steps, 5)])
     #plot(t2)
