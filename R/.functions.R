# Misc ####
agrestiCoull = function(X,n) {
  p=(X+1/2)/(n+1)
  return(sqrt(p*(1-p)/(n+1)))
}
erf = function(x) {
  2 * pnorm(x * sqrt(2)) - 1
}

printInt  = function(x1,x2, numDig=1)
  paste0('[',
         paste0(c(round(x1,numDig),round(x2,numDig)),
                collapse=','),
         ']')
estBS = function(error, nbs = 500, eps = 1) {

  mue = q95 = P1 = rmsd = mse = W = Wr = list()

  for (key in names(error)) {
    errors = error[[key]]

    # Mean Absolute Error
    bs = boot(errors,
              function(data, index) {
                mean(abs(data[index]))
              },
              nbs)
    mue[[key]] = prettyUnc(mean(bs$t), sd(bs$t), 1)

    # Mean Error
    bs = boot(errors,
              function(data, index) {
                mean(data[index])
              },
              nbs)
    mse[[key]] = prettyUnc(mean(bs$t), sd(bs$t), 1)

    # Root Mean Squared Deviation
    bs = boot(errors,
              function(data, index) {
                sd(data[index])
              },
              nbs)
    rmsd[[key]] = prettyUnc(mean(bs$t), sd(bs$t), 1)

    # 95th percentile
    bs = boot(errors,
              function(data, index) {
                quantile(abs(data[index]), probs = 0.95)
              },
              nbs)
    q95[[key]] = prettyUnc(mean(bs$t), sd(bs$t), 1)

    # P(Error < eps)
    bs = boot(errors,
              function(data, index) {
                sum(abs(data[index]) < eps) / length(index)
              },
              nbs)
    P1[[key]] = prettyUnc(mean(bs$t), sd(bs$t), 1)

    # Shapiro-Wilk normality statistics
    bs = boot(errors,
              function(data, index) {
                if (length(index) > 5000)
                  index = sample(index, 5000)
                shapiro.test(data[index])$statistic
              },
              nbs)
    W[[key]] = prettyUnc(mean(bs$t), sd(bs$t), 1)

    # Shapiro-Wilk normality p-value-based rejection
    bs = boot(errors,
              function(data, index) {
                if (length(index) > 5000)
                  index = sample(index, 5000)
                shapiro.test(data[index])$p
              },
              nbs)
    Wr[[key]] = signif(sum(bs$t < 0.01)/length(bs$t),2)
  }

  return(list(
    mue = mue,
    mse = mse,
    rmsd = rmsd,
    q95 = q95,
    P1 = P1,
    W = W,
    Wr = Wr
  ))

}
estBS0 = function(error, nbs = 500, eps = 1) {

  mue = q95 = P1 = rmsd = mse = W = Wr = list()
  umue = uq95 = uP1 = urmsd = umse = uW =  list()

  for (key in names(error)) {
    errors = error[[key]]

    # Mean Absolute Error
    bs = boot(errors,
              function(data, index) {
                mean(abs(data[index]))
              },
              nbs)
    mue[[key]] = mean(bs$t)
    umue[[key]] = sd(bs$t)

    # Mean Error
    bs = boot(errors,
              function(data, index) {
                mean(data[index])
              },
              nbs)
    mse[[key]] =mean(bs$t)
    umse[[key]] =sd(bs$t)

    # Root Mean Squared Deviation
    bs = boot(errors,
              function(data, index) {
                sd(data[index])
              },
              nbs)
    rmsd[[key]] = mean(bs$t)
    urmsd[[key]] = sd(bs$t)

    # 95th percentile
    bs = boot(errors,
              function(data, index) {
                quantile(abs(data[index]), probs = 0.95)
              },
              nbs)
    q95[[key]] = mean(bs$t)
    uq95[[key]] = sd(bs$t)

    # P(Error < eps)
    bs = boot(errors,
              function(data, index) {
                sum(abs(data[index]) < eps) / length(index)
              },
              nbs)
    P1[[key]] = mean(bs$t)
    uP1[[key]] = sd(bs$t)

    # Shapiro-Wilk normality statistics
    bs = boot(errors,
              function(data, index) {
                if (length(index) > 5000)
                  index = sample(index, 5000)
                shapiro.test(data[index])$statistic
              },
              nbs)
    W[[key]] = mean(bs$t)
    uW[[key]] = sd(bs$t)

    # Shapiro-Wilk normality p-value-based rejection
    bs = boot(errors,
              function(data, index) {
                if (length(index) > 5000)
                  index = sample(index, 5000)
                shapiro.test(data[index])$p
              },
              nbs)
    Wr[[key]] = signif(sum(bs$t < 0.01)/length(bs$t),2)
  }

  return(list(
    mue = mue,   umue = umue,
    mse = mse,   umse = umse,
    rmsd = rmsd, urmsd = urmsd,
    q95 = q95,   uq95 = uq95,
    P1 = P1,     uP1 = uP1,
    W = W,       uW = uW,
    Wr = Wr
  ))

}


dmue = function(X, index=1:nrow(X), uX = 0,...){
  v1 = abs(X[index,1])
  v2 = abs(X[index,2])
  N  = length(index)
  if(uX != 0) {
    pert = rnorm(N,0,uX) # Paired datasets
    v1 = v1 + pert
    v2 = v2 + pert
  }
  return(mean(v1)-mean(v2) )
}

estBS2 = function(error, uX = 0,
                  props = c('mue', 'wmue','mse', 'wmse',
                            'rmsd', 'q95hd', 'P1', 'W'),
                  est.corr = TRUE,
                  est.sip = TRUE,
                  nboot = 1000,
                  eps = 1) {
  options(boot.parallel = "multicore")

  if (class(error) == 'list') {
    n = names(error)
    error = as.matrix(as.data.frame(error))
    colnames(error) = n
  }

  wmse = function(X, index = 1:length(X), uX = 0) {
    x = X[index,]
    if(uX != 0) {
      ux2 = uX[index]^2
      t(
        apply(x,2,
              function(x) {
                sig2 = max( 0, var(x) - mean(ux2))
                w = 1 / (sig2 + ux2)
                Hmisc::wtd.mean(x,w)
              }
        )
      )
    } else {
      t(apply(x,2,mean,na.rm=TRUE))
    }
  }
  wmsece = function(X, index = 1:length(X), uX = 0) {
    x = X[index,]
    if(uX != 0) {
      ux = uX[index]
      w = 1/ux^2
      sw = sum(w)
      N = length(w)
      t(
        apply(x,2,
              function(x) {
                wm = weighted.mean(x,w)
                ud2ce = max(0,
                            (sum(w*(x-wm)^2) -(N-1)) /
                              (sw -sum(w^2)/sw)
                )
                wce = 1/(ud2ce+ux^2)
                wmsece= weighted.mean(x,wce)
              }
        )
      )
    } else {
      t(apply(x,2,mean,na.rm=TRUE))
    }
  }
  mse = function(X, index = 1:length(X), uX = 0) {
    x = X[index,]
    t(apply(x,2,mean,na.rm=TRUE))
  }
  rmsd = function(X, index = 1:length(X), uX = 0) {
    x = X[index,]
    if(uX != 0)
      x = x + rnorm(length(index),0,uX[index])
    t(apply(x,2,sd,na.rm=TRUE))
  }
  mue = function(X, index = 1:length(X), uX=0) {
    x = abs(X[index,])
    t(apply(x,2,mean,na.rm=TRUE))
  }
  wmue = function(X, index = 1:length(X), uX = 0) {
    x = X[index,]
    if(uX != 0) {
      ux2 = uX[index]^2
      t(
        apply(x,2,
              function(x) {
                sig2 = max( 0, var(x) - mean(ux2))
                w = 1 / (sig2 + ux2)
                Hmisc::wtd.mean(abs(x),w)
              }
        )
      )
    } else {
      x = abs(x)
      t(apply(x,2,mean,na.rm=TRUE))
    }
  }
  q95 = function(X, index = 1:length(X), uX = 0) {
    x = abs(X[index,])
    t(apply(x,2,quantile,probs=0.95,na.rm=TRUE))
  }
  wq95 = function(X, index = 1:length(X), uX = 0) {
    x = X[index,]
    if(uX != 0) {
      ux2 = uX[index]^2
      t(
        apply(x,2,
              function(x) {
                sig2 = max( 0, var(x) - mean(ux2))
                w = 1 / (sig2 + ux2)
                Hmisc::wtd.quantile(abs(x),w,0.95)
              }
        )
      )
    } else {
      x = abs(x)
      t(apply(x,2,quantile,probs=0.95,na.rm=TRUE))
    }
  }
  q95hd = function(X, index = 1:length(X), uX=0){
    x = abs(X[index,])
    t(apply(x,2,hd,q=0.95,na.rm=TRUE))
  }
  P1 = function(X, index = 1:length(X), uX = 0) {
    x = abs(X[index,])
    t(apply(x,2,function(x) sum(x < eps) / length(index)))
  }
  W = function(X, index = 1:length(X), uX = 0) {
    if (length(index) > 5000)
      index = sample(index, 5000)
    x = X[index,]
    t(apply(x,2,function(x) shapiro.test(x)$statistic))
  }

  # Process data
  results = list()
  methods = names(error)
  nm = length(methods)
  results[['props']] = props
  for (prop in props) {
    results[[prop]] = list()
    print(prop)
    statistic = get(prop) # associated function

    # Collectibe BS of all models to preserve correlations
    bs = boot::boot(error, statistic, R=nboot, uX = uX)

    # Unperturbed value for consistency check
    t0 = statistic(error,1:nrow(error),uX=0)

    # Store stats
    for (i in 1:nm) {
      m = methods[i]
      results[[prop]][['ref']][[m]] = t0[i]
      results[[prop]][['val']][[m]] = bs$t0[i]
      results[[prop]][['unc']][[m]] = sd(bs$t[,i])
      results[[prop]][['bs' ]][[m]] = bs$t[,i]
    }

    if(est.corr) {
      # Correlation of scores
      C = matrix(1, nrow = nm, ncol = nm)
      colnames(C) = rownames(C) = methods
      for (i in 1:(nm - 1)) {
        mi = methods[i]
        for (j in (i + 1):nm) {
          mj = methods[j]
          C[i, j] =
            cor(results[[prop]][['bs']][[mi]],
                results[[prop]][['bs']][[mj]],
                method = "spearman")
          C[j, i] = C[i, j]
        }
      }
      results[[prop]][['corr']] = C
    }
  }

  if(est.sip) {
    # Systematic improvement probability
    sip = usip = mg = umg = matrix(NA, nrow = nm, ncol = nm)
    for (i in 1:nm) {
      mi = methods[i]
      for (j in 1:nm) {
        if (j==i) next
        mj = methods[j]
        bs = boot::boot(
          cbind(error[[mi]], error[[mj]]),
          statistic = fsi,
          R = nboot,
          uX = uX)
        sip[i, j] = bs$t0[1]
        usip[i,j] = sd(bs$t[,1],na.rm=TRUE)
        mg[i, j]  = bs$t0[2]
        umg[i,j]  = sd(bs$t[,2],na.rm=TRUE)
      }
    }
    rownames(sip) =
      rownames(usip) =
      rownames(mg) =
      rownames(umg) =
      colnames(sip) =
      colnames(usip) =
      colnames(mg) =
      colnames(umg) = methods
    results[['sip']]  = sip
    results[['usip']] = usip
    results[['mg']]   = mg
    results[['umg']]  = umg
  }

  return(results)
}
disc = function(x,ux,y,uy,cxy=0) {
  abs(x-y)/sqrt(ux^2+uy^2-2*cxy*ux*uy)
}
cp = function(d,k=2) {
  ifelse(d<k,1,0)
}
cpInv = function(d,eps=0.05) {
  ifelse(d>eps,1,0)
}
calComp = function(S,prop,useCor=TRUE) {
  # Estimate compatibility matrix
  methList = names(S[[prop]]$val)
  nm = length(methList)
  D = matrix(0,nrow=nm,ncol=nm)
  for (i in 1:(nm-1)) {
    ni = methList[i]
    for(j in (i+1):nm) {
      nj = methList[j]
      x  = S[[prop]][['val']][ni]
      ux = S[[prop]][['unc']][ni]
      y  = S[[prop]][['val']][nj]
      uy = S[[prop]][['unc']][nj]
      cxy = 0
      if(useCor)
        cxy = S[[prop]][['corr']][ni,nj]
      D[i,j] = disc(x,ux,y,uy,cxy)
      D[j,i] = D[i,j]
    }
  }
  return(D)
}
calCompInv = function(S,prop) {
  # Estimate compatibility matrix based on inversions
  methList = names(S[[prop]]$val)
  nm = length(methList)
  D = matrix(NA,nrow=nm,ncol=nm)
  for (i in 1:nm) {
    ni = methList[i]
    for(j in 1:nm) {
      nj = methList[j]
      x  = S[[prop]][['bs']][[ni]]
      y  = S[[prop]][['bs']][[nj]]
      D[i,j] = sum(x<y)/length(x)
    }
  }
  return(D)
}

mcBS = function(N, score,
                nBS  = 1000,
                mean = c(0,0),
                sigma = diag(2),
                seed = NULL) {

  if(!is.null(seed))
    set.seed(seed)

  # BS estimate
  pg = del = rep(NA, nBS)
  for(irep in 1:nBS) {
    X = mvtnorm::rmvnorm(
      N,
      mean  = mean,
      sigma = sigma
    )
    bs = boot(
      X,
      statistic = fdif,
      R = 250,
      fscore = get(score)
    )
    del[irep] = bs$t0
    pg[irep]  = genpval(bs$t)
  }
  return(list(pg=pg,dif=del))
}
estPower = function(X1, X2, score = 'mse', nMC = 1000) {
# Estimate poxer of pg test usung
# normal approx. of datasets

  # Bivariate normal approx. of error sets
  p1 = appNorm(X1, score = score, nboot = nMC)
  p2 = appNorm(X2, score = score, nboot = nMC)
  rho = cor(X1, X2, method="spearman")
  means = c(p1$par[1], p2$par[1])
  sigma = matrix(c(p1$par[2]^2, rho*p1$par[2]*p2$par[2],
                   rho*p1$par[2]*p2$par[2], p2$par[2]^2),
                 nrow=2,ncol=2)
  print(means)
  print(sigma)
  N = length(X1)

  # Bootstrap within MC
  pg = mcBS(N, score,
            nBS   = nMC,
            mean  = means,
            sigma = sigma)

  mpg = mean(pg$pg)
  ppg = mean(pg$pg < 0.05)

  return(list(mpg = mpg, ppg = ppg, pgs = pg))

}




mse = function(X, index = 1:length(X)) {
  mean(X[index])
}
rmsd = function(X, index = 1:length(X)) {
  sd(X[index])
}
mue = function(X, index = 1:length(X)) {
  mean(abs(X[index]))
}
q95hd = function(X, index = 1:length(X)){
  # Quantile estimate by Harrell & Davis 1982
  hd(abs(X[index]), 0.95)
}
pinv = function (X,d0) {
  # Probability to have a sign different from d0's
  # The zeros (sign = 0) are excluded
  A = sum( sign(X) != sign(d0) )
  C = sum( X == 0 )
  (A - C)/length(X)
}
fcor = function(X, index=1:length(X),...){
  cor(X[index,1],X[index,2])
}
fcov = function(X, index=1:length(X),...){
  cov(X[index,1],X[index,2])
}
fdif = function(X, index=1:nrow(X),fscore,...){
  fscore(X[,1],index,...) - fscore(X[,2],index,...)
}
my5num = function(X) {
  c(
    hd(X, 0.05),
    hd(X, 0.25),
    hd(X, 0.5),
    hd(X, 0.75),
    hd(X, 0.95)
  )
}
muF = function(x) {
  mu=x[1]; sig=x[2]
  sig*sqrt(2/pi)*exp(-mu^2/(2*sig^2)) -mu*erf(-mu/(sqrt(2)*sig))
}
cdfF = function(x) {
  u = x[1]; mu = x[2]; sig = x[3]
  (erf((u+mu)/(sqrt(2)*sig))+erf((u-mu)/(sqrt(2)*sig)))/2
}
q95F_old = function(x,eps = seq(0,5,by=0.001)) {
  mu=x[1]; sig=x[2]
  G = apply(cbind(eps,mu,sig),1,cdfF)
  eps[(which(G>=0.95))[1]]
}
q95F = function(x) {
  mu=x[1]; sig=x[2]
  fz = function(x,mu,sig,prob) {
    cdfF(c(x,mu,sig)) - prob
  }
  mueF = muF(c(mu,sig))
  uniroot(f = fz, interval=c(mueF,mueF+6*sig),check.conv = TRUE,
          mu = mu, sig = sig, prob = 0.95)$root
}
fmod = function (pars,score) {
  mod = list()
  mod$mse  = mu  = pars[1]
  mod$rmsd = sig = pars[2]
  if( score == 'mue') {
    mod$score = muF(c(mu,sig))
  } else {
    mod$score = q95F(c(mu,sig))
  }
  return(mod)
}
fchi2 = function(pars,score,const,uconst) {
  # Fit the MSE and the target score
  mod = fmod(pars,score)
  chi2 =        (mod$mse   - const$mse)^2      #/ uconst$mse^2
  # chi2 = chi2 + (mod$rmsd  - const$rmsd)^2     #/ uconst$rmsd^2
  chi2 = chi2 + (mod$score - const[[score]])^2 #/ uconst[[score]]^2

  return(chi2)
}
appNorm = function(X, score = 'mue', nboot=1000) {
  # Find best normal approx with closest MSE, MUE, RMSD and Q95

  props = c('mse','rmsd',score)

  # Process data to get constraints
  const = uconst = list()
  for (prop in props) {
    bs = boot::boot(X, statistic = get(prop), R=nboot)
    const[[prop]]  = bs$t0
    uconst[[prop]] = sd(bs$t)
  }

  # Optimize
  startp = c(const$mse,const$rmsd)
  p = optim(
    par = startp,
    fn = fchi2,
    score = score,
    const = const,
    uconst = uconst)

  print(c('Const. :',signif(unlist(const),3)))
  print(c('Best   :',signif(unlist(fmod(p$par,score)),3)))

  return(p)
}
powPg = function (E1, E2, score = 'mue', nMC = 500, graph = FALSE) {
  if(graph) {
    par(mfrow=c(1,2))
    plot(ecdf(abs(E1)))
    plot(ecdf(abs(E2)), add=TRUE, col=2)
    abline(v = c(mue(E1),mue(E2)),col= 1:2)
    abline(v = c(q95hd(E1),q95hd(E2)),col= 1:2)
  }

  pg = estPower(E1, E2, score = score, nMC = nMC)

  if(graph) {
    hist(pg$pgs$pg,nclass=33)
    abline(v=median(pg$pgs$pg),col=2)
    abline(v=pg$mpg)
    abline(v=0.05,col=cols[3])
  }
  cat('<pg> = ', signif(pg$mpg,2),'\n')
  cat('P(pg<0.05) = ', signif(pg$ppg,2),'\n')
}



# Distribution of uncertainty on formation enthalpies ####
uncGen  = function(x) {
  rlnorm(1,
         meanlog = log(0.015 - 0.01*x),
         sdlog   = 0.45
  )
}
mcStats = function (X,Y,nMC = 100, numDig=1) {
  ## Random noise
  mue_t = q95_t = p1_t = matrix(NA, nMC ^ 2, 2)
  colnames(mue_t) = colnames(q95_t) = colnames(p1_t) = colnames(X)
  n = nrow(X)
  irun = 0
  for(iMC in 1:nMC) {
    uH = sapply(Y,uncGen)
    for(jMC in 1:nMC) {
      irun = irun +1
      for(key in colnames(X)) {
        AE = abs(rnorm(n,X[,key],uH))
        mue_t[irun,key] = mean(AE)
        q95_t[irun,key] = quantile(AE,probs = 0.95)
        p1_t[irun,key]  = sum(AE < 0.043)/length(AE)
      }
    }
  }
  ## Summary
  mue_m = apply(mue_t,2,function(x) mean(x,na.rm=TRUE))
  mue_s = apply(mue_t,2,function(x) sd(  x,na.rm=TRUE))
  q95_m = apply(q95_t,2,function(x) mean(x,na.rm=TRUE))
  q95_s = apply(q95_t,2,function(x) sd(  x,na.rm=TRUE))
  p1_m  = apply(p1_t ,2,function(x) mean(x,na.rm=TRUE))
  p1_s  = apply(p1_t ,2,function(x) sd(  x,na.rm=TRUE))

  # Format
  mue = q95 = p1 = list()
  for(key in colnames(X)) {
    mue[[key]] = prettyUnc(mue_m[[key]],mue_s[[key]], numDig=numDig)
    q95[[key]] = prettyUnc(q95_m[[key]],q95_s[[key]], numDig=numDig)
    p1[[key]]  = prettyUnc(p1_m[[key]], p1_s[[key]] , numDig=numDig)
  }
  return(list(mue=mue, q95=q95, p1=p1))
}

# Iterative Re-Weighted Least Squares ####

mpuReg <- function(x, y, uy = NULL, Np = 2) {
  # Estimate mean prediction uncertainty
  # after linear regression by IRWLS
  #
  # Linear fit by a model : yi = f(xi) + d + ei
  # with ei ~ N(0,uyi) : experimental uncertainty
  # and  d ~ N(0,ud)   : model error uncertainty

  if(is.null(uy))
    stop('>>> mpuReg: uy should be a numeric (>0) vector <<<')
  uy2 = uy^2

  # Build Design matrix
  X = rep(1, length(x))
  for (i in 2:Np)
    X = cbind(X, x ^ (i - 1))

  # First pass to estimate ud2
  S = diag(uy2, nrow = length(uy2))
  Sm1 = solve(S)
  beta0 = solve(t(X) %*% Sm1 %*% X) %*% (t(X) %*% Sm1 %*% y)
  bm = X %*% beta0
  # Compute unexplained variance in residuals
  ud2 = max( 0, var(y - bm) - mean(uy2))

  # Second pass with new weights
  S = diag(uy2 + ud2, nrow = length(uy2))
  Sm1 = solve(S)
  beta0 = solve(t(X) %*% Sm1 %*% X) %*% (t(X) %*% Sm1 %*% y)
  bm = X %*% beta0
  ud2 = max( 0, var(y - bm) - mean(uy2))

  # Additional passes to converge reweighting scheme
  # Useless if uy is homogeneous
  if(sd(uy) != 0) {
    for (i in 1:5) {
      # Redefine weights after variance partition
      S = diag(uy2  + ud2, nrow = length(uy2))
      Sm1 = solve(S)
      beta0 = solve(t(X) %*% Sm1 %*% X) %*% (t(X) %*% Sm1 %*% y)
      bm = X %*% beta0
      ud2 = max( 0, var(y - bm) - mean(uy2))
    }
  }

  # Params covariance / Laplace approx.
  Sbeta = solve(t(X) %*% Sm1 %*% X)

  # Covariance of predictions (confidence)
  SP = X %*% Sbeta %*% t(X)

  uP2 = diag(SP) + ud2

  # Mean prediction uncertainty
  mpu = sqrt(mean(uP2))

  return(
    list(
      coefficients  = beta0,
      coef.vcov     = Sbeta,
      fitted.values = bm,
      fitted.unc    = sqrt(uP2),
      d             = sqrt(ud2),
      mpu           = mpu
    )
  )
}
betaf = function(M,
                 index = 1:nrow(M),
                 Np ) {
  x  = M[index, 1]
  y  = M[index, 2]
  uy = M[index, 3]
  reg = mpuReg(x, y, uy, Np = Np)
  reg$coefficients
}
uPf = function(M,
               index = 1:nrow(M),
               Np ) {
  x  = M[index, 1]
  y  = M[index, 2]
  uy = M[index, 3]
  reg = mpuReg(x, y, uy, Np = Np)
  io = order(x)
  c(x[io],reg$fitted.unc[io])
}
mpuf = function(M,
                index = 1:nrow(M),
                Np ) {
  x  = M[index, 1]
  y  = M[index, 2]
  uy = M[index, 3]
  reg = mpuReg(x, y, uy, Np = Np)
  reg$mpu
}
irw_reglin = function(x,y,uy){
  n=length(x)

  # Define weights
  w=1/uy^2

  d = sum(w)*sum(w*x^2)-sum(w*x)^2
  b = (sum(w)*sum(w*x*y)-sum(w*x)*sum(w*y))/d
  a = sum(w*y)/sum(w)-b*sum(w*x)/sum(w)
  x2= sum(w*(y-a-b*x)^2)
  if(x2 <= qchisq(0.05,df=n-2))
    stop(paste('exp. errors too large',x2))

  for (i in 1:5) {

    ud2 = max( 0,
               1 /(n-2) *sum((y-a-b*x)^2)-mean(uy^2)
    )

    # Redefine weights after variance partition
    w=1/(ud2+uy^2)

    d = sum(w)*sum(w*x^2)-sum(w*x)^2
    b = (sum(w)*sum(w*x*y)-sum(w*x)*sum(w*y))/d
    a = sum(w*y)/sum(w)-b*sum(w*x)/sum(w)
    ub2 =  sum(w)/d
    ua2 =  sum(w*x^2)/d
    cab = -sum(w*x)/sum(w)*ub2
  }

  return(list(a=a,b=b,ud2=ud2,ua2=ua2,ub2=ub2,cab=cab))
}
irw_predlin = function(r,x) {
  a = r$a; b=r$b; ud2= r$ud2; ua2=r$ua2; ub2=r$ub2; cab=r$cab
  p   = a + b*x
  ul2 = ua2 + x^2*ub2 + 2*x*cab
  up2 = ud2 + ul2
  return(list(p=p,up=up2^0.5,ul=ul2^0.5,ratio=ul2/up2))
}

irw_lm = function(x, y, uy = NULL, mod = 'y~1+x') {
  # Linear fit by a model
  # yi = f(xi) + d + ei
  # with ei ~ N(0,uyi) : experimental uncertainty
  # and  d ~ N(0,ud)   : model error uncertainty

  if(is.null(uy))
    stop('>>> irw_lm: uy should be a numeric (>0) vector <<<')

  # Regression
  reg = lm(mod,
           data = data.frame(x=x,y=y,uy=uy),
           weights = 1/uy^2)
  res = reg$residuals

  # Check reduced chi-squared
  chi2 = sum((res/uy)^2)
  nu   = reg$df.residual
  if (chi2 <= qchisq(0.05, df = nu))
    stop(paste0('>>> irw_lm: exp. errors too large; chi2 = ',chi2,' <<<'))

  # Compute unexplained variance in residuals
  ud2 = max( 0, var(res) - mean(uy^2))

  # Loop to converge reweighting scheme
  # Useless if uy is homogeneous
  if(sd(uy) != 0) {
    for (i in 1:5) {
      # Redefine weights after variance partition
      reg = lm(mod,
               data = data.frame(x=x,y=y,uy=uy),
               weights = 1/(ud2 + uy ^ 2))
      res = reg$residuals
      ud2 = max( 0, var(res) - mean(uy^2))
    }
  }
  return(
    list(
      reg = reg,
      ud2 = ud2,
      uy  = uy
    )
  )
}
irw_predlin_lm = function(r,x) {

  reg = r$reg
  ud2 = r$ud2
  uy  = r$uy

  pr = predict.lm(
    reg,
    newdata=data.frame(x=x),
    se.fit = TRUE,
    weights = reg$weights,
    scale = sqrt(mean(uy^2))
  )
  p   = pr$fit
  up2 = ud2 + pr$se.fit^2
  return(list(p=p,up=up2^0.5,ud=ud2^0.5,ratio=ud2/up2))
}

irw_ls_lm=function(x,y,uy=NULL, mod = 'y~1+x', fac=2) {

  r = irw_lm(x,y,uy,mod)

  # Mean variance and max contrib of line uncert.
  ns = 1000
  x0 = seq(min(x),max(x),len=ns)
  pr  = irw_predlin_lm(r,x0)
  vmean = mean(pr$up^2)
  ratio = max(pr$ratio)

  # Coverage score
  pr = irw_predlin_lm(r,x)
  vtot=mean(uy^2)+vmean
  score = sum( y>=(pr$p-fac*vtot^0.5) & y<=(pr$p+fac*vtot^0.5) )
  score = score/length(x)

  # LOO stat
  yloo=c()
  for (iloo in 1:length(y)) {
    y0 = y[-iloo]
    x0 = x[-iloo]
    uy0= uy[-iloo]
    r0 = irw_lm(x=x0,y=y0,uy=uy0,mod)
    yloo[iloo] = irw_predlin_lm(r0,x[iloo])$p
  }
  PSS = sum((y - yloo   )^2)
  TSS = sum((y - mean(y))^2)
  Q2=1 - PSS/TSS

  return(list(reg = r, ud=r$ud2^0.5,
              vmean=vmean, ratioMax=ratio,
              score=score, Q2=Q2))
}

# N=100
# eps =4
# x = 1:N
# xp = seq(0,2*max(N),by=0.1)
# mod = 'y~1+x'
# y = eps*rnorm(N)
#
# par(mfrow=c(1,2))
# uy = rep(eps,N)
# C = irw_lm(x,y,uy,mod)
# print( 1/C$reg$df.residual * sum(C$reg$residuals^2/(uy^2+C$ud2)))
# P = irw_predlin_lm(C,xp)
#
# plot(x,y, main =mean(P$ratio),ylim = 4*eps*c(-1,1))
# grid()
# segments(x,y-2*uy,x,y+2*uy)
# lines(xp,P$p,col=4)
# lines(xp,P$p-2*P$up,col=4)
# lines(xp,P$p+2*P$up,col=4)
# lines(xp,P$p-2*sqrt(P$up^2+mean(uy^2)),col=4,lty=2)
# lines(xp,P$p+2*sqrt(P$up^2+mean(uy^2)),col=4,lty=2)
#
# uy = rep(eps/10,N)
# C = irw_lm(x,y,uy,mod)
# print( 1/C$reg$df.residual * sum(C$reg$residuals^2/(uy^2+C$ud2)))
# P = irw_predlin_lm(C,xp)
# plot(x,y, main =mean(P$ratio),ylim = 4*eps*c(-1,1))
# grid()
# segments(x,y-2*uy,x,y+2*uy)
# lines(xp,P$p,col=4)
# lines(xp,P$p-2*P$up,col=4)
# lines(xp,P$p+2*P$up,col=4)
# lines(xp,P$p-2*sqrt(P$up^2+mean(uy^2)),col=4,lty=2)
# lines(xp,P$p+2*sqrt(P$up^2+mean(uy^2)),col=4,lty=2)




irw_ls=function(x,y,uy,fac=1) {

  r = irw_reglin(x,y,uy)

  # Mean variance and max contrib of line uncert.
  ns = 1000
  x0 = seq(min(x),max(x),len=ns)
  pr  = irw_predlin(r,x0)
  vmean = mean(pr$up^2)
  ratio = max(pr$ratio)

  # Coverage score
  pr = irw_predlin(r,x)
  vtot=mean(uy^2)+vmean
  score = sum( y>=(pr$p-fac*vtot^0.5) & y<=(pr$p+fac*vtot^0.5) )
  score = score/length(x)

  # LOO stat
  yloo=c()
  for (iloo in 1:length(y)) {
    y0 = y[-iloo]
    x0 = x[-iloo]
    uy0= uy[-iloo]
    r0 = irw_reglin(x0,y0,uy0)
    yloo[iloo] = irw_predlin(r0,x[iloo])$p
  }
  PSS = sum((y - yloo   )^2)
  TSS = sum((y - mean(y))^2)
  Q2=1 - PSS/TSS

  return(list(a=r$a, b=r$b, ud=r$ud2^0.5,
              vmean=vmean, ratioMax=ratio,
              score=score, Q2=Q2,
              ua2=r$ua2, ub2=r$ub2, cab=r$cab))
}
# Chemical Composition analysis ####
getElts = function(sys) {
  sys1 = sub('Oa','O',sys)
  m = gregexpr('([[:upper:]][[:lower:]]?)',sys1)
  elements = sort(unique(unlist(regmatches(sys1,m))))
  return(elements)
}
getComposBi = function(sys) {
  # Get elements in binary systems
  sys1 = sub('Oa','O',sys)
  m = gregexpr('([[:upper:]][[:lower:]]?)',sys1)
  compos = regmatches(sys1,m)
  return(matrix(unlist(compos),ncol=2,byrow=TRUE))
}
# Plot ####
plotXY = function(X, Y, dY, gPars, units = 'meV'){
  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(mfrow=c(1,1),mar=mar,mgp=mgp,
      pty=pty,tcl=tcl, cex=cex, lwd=lwd)


  xlim = ylim = range(c(X,Y))
  plot(X, Y, pch=16, cex = 0.5, col=cols_tr2[5], log='',
       xlab = paste0('Calc. Enthalpy ',units),
       xlim = xlim,
       ylab = paste0('Experiment ',units),
       ylim = ylim,
       main='')
  grid()
  abline(a=0, b = 1, lty=3)
  abline(lm(Y~X),col=cols[5])
  box()
}







plotZscoreQqnorm <- function(R,
                             sig,
                             lim = NULL,
                             title = '',
                             gPars) {

  # Monte Carlo estimation of 95% CI of qqnorm
  # for a given sample size N
  N   = length(R)
  nMC = 1000
  uly = matrix(0, ncol = N, nrow = nMC)
  for (i in 1:nMC) {
    t = sort(rnorm(N))
    qt = qqnorm(t, plot.it = FALSE)
    uly[i, ] = qt$y
  }
  q95 = t(
    apply(
      uly, 2,
      quantile,
      names = FALSE,
      probs = c(0.025, 0.975)
    )
  )
  rm(uly)
  ulx = sort(qt$x) # Set of quantiles for N points

  # Z-scores
  Zs = R / sig

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend = 2
  )

  if (is.null(lim))
    lim = max(3, max(abs(Zs)))

  q = qqnorm(
    Zs,
    main = '',
    pch = 16,
    col = cols[5],
    xlim = lim * c(-1, 1),
    ylim = lim * c(-1, 1)
  )
  grid()
  polygon(c(ulx, rev(ulx)),
          c(q95[, 1], rev(q95[, 2])),
          border = NA,
          col = cols_tr2[3])
  points(q$x,
         q$y,
         pch = 16,
         cex = 0.8,
         col = cols[5])
  qqline(Zs, col = 2)
  abline(a = 0, b = 1)
  out = (Zs - q95[, 1]) * (Zs - q95[, 2]) > 0
  text(x = ulx[out],
       y = Zs[out],
       labels = names(R)[out])
  legend(
    'topleft',
    inset = 0.025,
    title = title,
    title.adj = 0,
    bty = 'n',
    legend = c('Data',
               '95% CI'),
    pch = c(16, -1),
    col = c(cols[5], cols_tr2[3]),
    lty = c(NA, 1),
    lwd = c(0, 30)
  )
  box()
}

plotMcECDF = function (X,
                   Y,
                   xlab = '|Error| (eV/atom)',
                   xmax = NULL,
                   title = '',
                   show.leg = TRUE,
                   col.index = 1:ncol(X),
                   nMC = 100,
                   gPars) {

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend= 2
  )

  ## Random noise
  irun = 0
  n = nrow(X)
  prob = (1:n)/n
  tab = array(NA,dim = c(n,nMC^2,2),
              dimnames = list(NULL,NULL,colnames(X)) )
  for(iMC in 1:nMC) {
    uH = sapply(Y,uncGen)
    for(jMC in 1:nMC) {
      irun = irun +1
      for(key in colnames(X)) {
        samp = abs(rnorm(n,X[,key],uH))
        tab[,irun,key] = sort(samp)
      }
    }
  }
  stat = array(NA,dim = c(n,3,2),
               dimnames = list(NULL,NULL,colnames(X)) )
  for(key in colnames(X)) {
    stat[,,key] = t(
      apply(
        tab[,,key],1,
        quantile,
        names = FALSE,
        probs=c(0.025,0.5,0.975)
      )
    )
  }

  for (icol in 1:ncol(X)) {
    if (icol == 1) {
      plot(
        stat[,2,icol],
        prob,
        type = 'l',
        col  = cols[col.index[icol]],
        xlab = xlab,
        xlim = c(0, xmax),
        main = '',
        yaxs = 'i',
        ylab = 'Probability'
      )
      grid(lwd = 2)
    } else {
      lines(stat[,2,icol], prob, col = cols[col.index[icol]])
    }

    # sigp = sqrt(prob * (1 - prob) / length(x))
    polygon(c(stat[,1,icol], rev(stat[,3,icol])),
            c(prob, rev(prob)),
            col = cols_tr2[col.index[icol]],
            border = NA)
    iq = which(prob >= 0.95)[1]
    abline(v = stat[iq, 2, icol],
           col = cols[col.index[icol]],
           lty = 2)
  }

  abline(h=0.95,col=2,lty=2)
  mtext(text = '0.95 ', at=0.95, side = 2,
        col=2, cex=cex, las=2)
  box()
  if(show.leg)
    legend('bottomright', title = title, title.adj = 0,
           legend = colnames(X), bty='n',
           col= cols_tr2[col.index], lty=1, lwd =30,
           cex = cex.leg)

}

plotMcDeltaECDF <- function(X, Y,
                            meth1,
                            meth2,
                            eps = NULL,
                            xmax = NULL,
                            xlab = 'Delta Abs. errors',
                            nMC = 50,
                            gPars) {
  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend = 2,
    xaxs = 'i',
    yaxs = 'i'
  )

  ## Random noise
  irun = 0
  n = nrow(X)
  prob = (1:n)/n

  tab = array(NA,dim = c(n,nMC^2,1))
  for(iMC in 1:nMC) {
    uH = sapply(Y,uncGen)
    for(jMC in 1:nMC) {
      irun = irun +1
      pert = (rnorm(n,0,uH))
      samp = abs(X[,meth1]+pert) - abs(X[,meth2]+pert)
      tab[,irun,1] = sort(samp)
    }
  }
  stat = array(NA,dim = c(n,3,1))
  stat[,,1] = t(
    apply(
      tab[,,1],1,
      quantile,
      names=FALSE,
      probs=c(0.025,0.5,0.975)
    )
  )

  if(is.null(xmax))
    xmax = max(X)

  plot(
    stat[,2,1],
    prob,
    type = 'l',
    col  = cols[5],
    xlab = xlab,
    # xlim = c(0, xmax),
    main = '',
    yaxs = 'i',
    ylab = 'Probability'
  )
  grid(lwd = 2)
  polygon(c(stat[,1,1], rev(stat[,3,1])),
          c(prob, rev(prob)),
          col = cols_tr2[5],
          border = NA)

  if (!is.null(eps))
    polygon(
      x = c(-eps, -eps, eps, eps),
      y = c(0, 1, 1, 0),
      col = cols_tr2[3],
      border = NA
    )
  abline(v = 0,lty=2)
  box()
}
plotOnPeriodicTable = function(X,info,label='') {
  # PLot values of X on periodic table
  # X should be a named vector with names as periodic table symbols

  # Script from https://chepec.se/2014/11/16/element-data/
  # Ref. Ahmed, Taha. 2014. “Properties of the elements.” https://doi.org/10.6084/m9.figshare.1295585
  #library(ggplot2)
  #library(dplyr)
  #library(grid)

  # load(url("http://files.figshare.com/1873544/element_data.rda"))
  load("../data/element_data.rda") # local copy


  b = rep(NA,length(values$Symbol))
  names(b) = values$Symbol
  b[names(X)] = X
  values$b = b

  labs = rep(NA,length(values$Symbol))
  names(labs) = values$Symbol
  labs[names(info)] = info

  ticks = pretty(b)

  pl = ggplot() +
    geom_point(data = values,
               # size 14 fits well with fig.width = 9, fig.height = 5.25
               # size = 14,
               size = 8,
               # shape #22 allows both fill and colour to be
               # shape #15 only registers colour (is always filled)
               # shape #0 only registers colour (is always empty)
               shape = 15,
               aes(y = Graph.Period, x = Graph.Group, colour = b)) +
    geom_text(data = values,
              colour = "white",
              # size = 3,
              size=2,
              fontface = "bold",
              aes(label = Symbol, y = Graph.Period-0.25, x = Graph.Group-0.25)) +
    geom_text(data = values,
              colour = "white",
              # size = 4,
              size = 2,
              fontface = "bold",
              aes(label = labs, y = Graph.Period+0.17, x = Graph.Group)) +
    scale_x_continuous(breaks = seq(min(values$Graph.Group),
                                    max(values$Graph.Group)),
                       limits = c(min(values$Graph.Group) - 1,
                                  max(values$Graph.Group) + 1),
                       expand = c(0,0)) +
    scale_y_continuous(trans = "reverse",
                       breaks = seq(min(values$Graph.Period),
                                    max(values$Graph.Period)),
                       limits = c(max(values$Graph.Period) + 1,
                                  min(values$Graph.Period) - 1.5),
                       expand = c(0,0)) +
    # set breaks and labels in the colourbar legend
    scale_colour_gradient2(breaks = ticks,
                           labels = ticks,
                           limits = c(0,1),
                           # low = muted('blue'),
                           mid = 'gray70',
                           # high = muted('red'),
                           na.value = "grey90") +
    # plot title (usually property and unit)
    # annotate("text", x = 8, y = 0.6,
    #          vjust = 0,
    #          label = label,
    #          # label = paste("Density/", "10^3*~kg~m", "^-3", sep=""),
    #          # parse required if using superscripts or subscripts
    #          parse = FALSE) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # plot.margin = unit(c(0, 0, -0.85, -0.85), "line"),
          plot.margin = unit(-0.5*c(1, 1, 1, 1), "line"),
          aspect.ratio = 5.25/9,
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          # center (roughly) over transition metal block
          legend.position = c(0.42, 0.91),
          legend.justification = c(0.5, 1),
          legend.direction = "horizontal",
          legend.key.height = unit(0.4, "line"),
          # make the legend colourbar a little longer
          # legend.key.width = unit(2.5, "line"),
          legend.key.width = unit(1.5, "line"),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"))

  return(pl)
}



myVioplot = function (datas, range = 1.5, h = NULL,
                      ylim = NULL, names = NULL,
                      horizontal = FALSE,
                      col = 1:length(datas),
                      border = "black", lty = 1,
                      lwd = 1, rectCol = cols_tr[5],
                      colMed = "white", pchMed = 19,
                      at, add = FALSE, wex = 1,
                      drawRect = TRUE, side='b',...)
{
  #   datas <- list(x, ...)
  n <- length(datas)
  if (missing(at))
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h)))
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.975, na.rm=TRUE)
    q3[i] <- quantile(data, 0.025, na.rm=TRUE)
    med[i] <- median(data, na.rm=TRUE)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i],
                                               data.max))
    #     smout <- do.call("sm.density", c(list(data, xlim = est.xlim),
    #                                      args))
    #     hscale <- 0.4/max(smout$estimate) * wex
    #     base[[i]] <- smout$eval.points
    #     height[[i]] <- smout$estimate * hscale
    smout <- do.call("density", args=list(data,cut=0))
    hscale <- 0.4/max(smout$y) * wex
    base[[i]] <- smout$x
    height[[i]] <- smout$y * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1)
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.5 * wex
  if (!add)
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim, ...)
      #       axis(2)
      ###       axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      if( side =='b' ) {
        x = c(at[i] - height[[i]], rev(at[i] + height[[i]]))
        y = c(base[[i]], rev(base[[i]]))
      } else {
        if(side== 'l') {
          x = c(at[i] - height[[i]], rev(at[i]+ 0*height[[i]]))
          y = c(base[[i]], rev(base[[i]]))
        } else {
          x = c(at[i]- 0*height[[i]], rev(at[i] + height[[i]]))
          y = c(base[[i]], rev(base[[i]]))
        }
      }
      polygon(x,y, col = col[i], border = border,
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd,
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2,
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim,...)
      #       axis(1)
      ###       axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {   if( side =='b' ) {
      x = c(at[i] - height[[i]], rev(at[i] + height[[i]]))
      y = c(base[[i]], rev(base[[i]]))
    } else {
      if(side== 'l') {
        x = c(at[i] - height[[i]], rev(at[i]+ 0*height[[i]]))
        y = c(base[[i]], rev(base[[i]]))
      } else {
        x = c(at[i]- 0*height[[i]], rev(at[i] + height[[i]]))
        y = c(base[[i]], rev(base[[i]]))
      }
    }
      polygon(y,x,
              col = col[i], border = border,
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd,
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] +
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med,
                 q1 = q1, q3 = q3))
}
plot2Dens = function(data1,data2,
                     col1='green',
                     legends=NULL,
                     x.leg='topleft',
                     y.leg=NULL) {
  ylim=range(data1)
  ll=list()
  for (i in 1:ncol(data1))
      ll[[i]] = data1[,i]
  myVioplot(ll,col=rep(cols[5],length(ll)),
            drawRect=TRUE, horiz=FALSE,
            yaxt='n',ylab='',ylim=ylim,
            xaxt='n',xlab='',
            side='l',lwd=2)
  grid()
  # axis(1,lwd=4,padj=-0.7)
  # axis(3,lwd=4,padj= 0.7)

  ll=list()
  for (i in 1:ncol(data2))
    ll[[i]] = data2[,i]
  myVioplot(ll,col=rep(cols[2],length(ll)),
            drawRect=TRUE, horiz=FALSE,
            yaxt='n',ylab='',ylim=ylim,
            xaxt='n',xlab='',
            side='r',lwd=2,
            add=TRUE)
  # title(xlab= paste0('Errors on ',
  #                    tolower(casesLong[ic]),
  #                    casesUnits[ic]),line=1.5)
  # mtext(rev(paste0(names,' ')),side=2,
  #       at=(1:length(names))-0.5/heightFactor,las=2,
  #       srt=90,adj=1,padj=0,cex=8,col='black')


  # grid(ny=0,col='gray30',lwd=2,lty=3)
  # abline(v=0,lwd=2)
  # box(lwd=6)
  # if(!is.null(legends))
  #   legend(x=x.leg, y=y.leg, cex=0.8,
  #          legend = legends, bty='y', bg='gray90',
  #          density = -1, fill=c(col1,col2), text.col=1)
}

# Unchecked ####

plotErrors = function(X, Y, dY, gPars, units = 'meV'){
  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(mfrow=c(1,1),mar=mar,mgp=mgp,
      pty=pty,tcl=tcl, cex=cex, lwd=lwd)

  err = X - Y
  ylim=range(c(err-2*dY,err+2*dY))
  plot(X, err, pch=16, col=cols[2], log='',
       xlab = paste0('Calculated value / ',units),
       ylab = paste0('Error / ',units),
       ylim = ylim,
       main='')
  grid()
  segments(X,err-2*dY,X,err+2*dY,col=cols[2])
  abline(h=0, col = 1)
  abline(lm(err~0+X),lty=2,col=4)
  box()
}
plotLocCI <- function(X,R,uy,nSeg=10, gPars, ylim=NULL) {
  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))
  par(mfrow=c(1,1),mar=mar,mgp=mgp,
      pty=pty,tcl=tcl, cex=cex, lwd=lwd)

  # del =  diff(range(X)) / nSeg # Width of segment


  # Take packets of equal size (afap)
  io = order(X)
  X  = X[io]; R = R[io]; uy = uy[io]
  np = length(R)/nSeg

  lci  = -1.96 * uy
  uci  =  1.96 * uy

  prop = uprop = c()
  for (j in 1:nSeg) {
    # sel =  (X0 + del * (j - 1) <= X) & (X < X0 + del * j)
    sel =  ((j - 1) * np + 1):(j * np)
    if (j == nSeg)
      sel =  ((j - 1) * np + 1):length(R)
    c   = sum((lci[sel] - R[sel]) * (uci[sel] - R[sel]) < 0)
    p = c / length(sel)
    prop[j]  = p
    uprop[j] = (p * (1 - p) / length(sel)) ^ 0.5
  }

  # Coverage statistics
  if(is.null(ylim))
    ylim = c(min(prop),1.0)
  plot(1:nSeg,prop, ylim=ylim, col=cols[5], pch=16,
       xlab = 'X-segment #',
       ylab = 'Probability',
       main = 'Coverage probability of 95% PUI-CI')
  grid()
  segments(
    1:nSeg,pmax(0,prop-1.96*uprop),
    1:nSeg,pmin(1,prop+1.96*uprop),
    col = cols[5])
  abline(h=0.95,col=cols[3],lty=2)
}
plotZscoreEcdf <- function(R, sig, gPars) {
  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  # Estimate absolute Zscore
  N      = length(R)
  Zscore = R/sig
  Z      = abs(Zscore)

  # ECDF of abs. Zscore
  xZcdf  = sort(Z)
  yZcdf  = seq(0,1,length.out = N)

  # Diff. of ECDF with CDF of Half-Normal dist.
  Zdif   = yZcdf - erf(xZcdf/sqrt(2))

  # CDF of Half-Normal dist. and approx. sampling uncertainty
  X     = seq(0,max(xZcdf),by=0.001)
  Ncdf  = erf(X/sqrt(2))
  uNcdf = sqrt(Ncdf*(1-Ncdf)/N)

  # Plot
  par(mfrow=c(1,1),mar=mar,mgp=mgp,
      pty=pty,tcl=tcl, cex=cex, lwd=lwd)

  plot(xZcdf, erf(xZcdf/sqrt(2)), type ='l', col=cols[2],
       xlab = expression("|"~Z[score]~"|"),
       ylab = 'CDF',
       ylim = c(0,1), yaxs = 'i')
  grid(lwd=1)
  polygon(
    c(X,rev(X)),
    c(Ncdf + 3 *uNcdf,rev(Ncdf - 3 *uNcdf)),
    col=cols_tr2[5], border=NA
  )
  points(xZcdf, yZcdf, col=cols[2],pch=16)
  box()

}
plotZscoreOnion <- function(R, sig, gPars) {
  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  # Estimate absolute Zscore
  N      = length(R)
  Zscore = R/sig
  Z      = abs(Zscore)

  # ECDF of abs. Zscore
  xZcdf  = sort(Z)
  yZcdf  = seq(0,1,length.out = N)

  # Diff. of ECDF with CDF of Half-Normal dist.
  Zdif   = yZcdf - erf(xZcdf/sqrt(2))

  # CDF of Half-Normal dist. and approx. sampling uncertainty
  X     = seq(0,max(xZcdf),by=0.001)
  Ncdf  = erf(X/sqrt(2))
  uNcdf = sqrt(Ncdf*(1-Ncdf)/N)

  #Plot
  par(mfrow=c(1,1),mar=mar,mgp=mgp,
      pty=pty,tcl=tcl, cex=cex, lwd=lwd)

  ylim = 1.1*max(c(3*uNcdf,abs(Zdif)))
  plot(xZcdf,Zdif,type='n',
       xlab = expression("|"~Z[score]~"|"),
       ylab = expression(Delta[CDF]),
       ylim = ylim*c(-1,1))
  grid(lwd=1)

  # Contour of 99% CI
  polygon(c(X,rev(X)),3*c(uNcdf,-rev(uNcdf)),col=0)
  # Fill-in with color
  for(i in seq(0.5,3,by=0.5))
    polygon(c(X,rev(X)),i*c(uNcdf,-rev(uNcdf)),
            col=cols_tr[5],border=NA)

  # Data points
  points(xZcdf,Zdif,col=cols[2],pch=16)

  # Colorize outliers
  t = erf(xZcdf/sqrt(2))
  over = Zdif > 3*sqrt(t*(1-t)/N)
  points(xZcdf[over],Zdif[over],col=cols[1],pch=16)
  below = Zdif < -3*sqrt(t*(1-t)/N)
  points(xZcdf[below],Zdif[below],col=cols[1],pch=16)

  # Percentage of outliers
  outPerc = 100 * sum(over | below) / N
  legend('topright',
         title=paste0('Outliers: ',signif(outPerc,2),'%'),
         legend=NA, bty='n')
  box()
}
# Misc. functions
resid = function(p, D, X) {
  D - model(p, X)
}
nLogLik = function(p, D, X) {
  0.5 * sum(resid(p, D, X) ^ 2 / uRef ^ 2)
}
covProb <- function(X,
                    D,
                    uD,
                    p0,
                    nMC = 500,
                    nMCUP = 1000,
                    nSeg = 10,
                    uRef = 0) {
  del =  diff(range(X)) / nSeg # Width of segment

  tab = matrix(NA, nrow = nMC, ncol = nSeg)
  for (imc in 1:nMC) {
    D   = model(p0, X) + rnorm(length(X), 0, uD)
    reg = optim(
      par = p0,
      fn = nLogLik,
      D = D,
      X = X,
      hessian = TRUE
    )
    opt = reg$par
    sig = solve(reg$hessian)

    # Variance of residuals
    R   = resid(opt, D, X)
    v1  = var(R)

    # MC-UP: var of model
    S_in  = mvtnorm::rmvnorm(nMCUP, opt, sig)
    S_out = t(apply(S_in, 1, function(p)
      model(p, X)))
    vp    = apply(S_out, 2, var)
    v2    = mean(vp)

    # Model confidence interval (95%)
    up0   = sqrt(vp)
    lci0  = -1.96 * up0
    uci0  =  1.96 * up0

    # Model prediction interval (95%)
    lci1  = -1.96 * sqrt(up0^2 + uRef^2)
    uci1  =  1.96 * sqrt(up0^2 + uRef^2)

    # Inflation factor
    puiF = sqrt(v1 / v2)

    # PUI confidence interval (95%)
    up   = up0 * puiF
    lci  = -1.96 * up
    uci  =  1.96 * up

    # Estimate coverage probability
    prop = c()
    X0 = min(X)
    for (j in 1:nSeg) {
      sel     =  (X0 + del * (j - 1) <= X) & (X < X0 + del * j)
      c       = sum((lci[sel] - R[sel]) * (uci[sel] - R[sel]) < 0)
      prop[j] = c / sum(sel)
    }
    tab[imc, ] = prop
  }
  m   = colMeans(tab)
  q95 = t(
    apply(
      tab, 2,
      quantile,
      names = FALSE,
      probs = c(0.025, 0.975)
    )
  )

  return(list(
    m = m,
    q95 = q95,
    lci0 = lci0,
    uci0 = uci0,
    lci1 = lci1,
    uci1 = uci1,
    lci = lci,
    uci = uci
  ))
}
plotCovProb <- function(X,R,uy,cp,gPars) {
  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend= 2
  )


  ylim = range(R+2*uy, R-2*uy, cp$lci, cp$uci)
  plot( X, R,
        type ='n',
        ylim = ylim,
        pty = 's',
        xaxs = 'i',
        yaxs = 'i',
        xlab = 'X',
        ylab = 'Residuals',
        main = ''
  )
  grid()
  abline(h = 0)
  segments(X,R-2*uy,X,R+2*uy,col = cols[5])
  polygon(c(X, rev(X)),
          c(cp$lci0, rev(cp$uci0)),
          col = cols_tr2[2], border = NA)
  polygon(c(X, rev(X)),
          c(cp$lci1, rev(cp$uci1)),
          col = cols_tr2[3], border = NA)
  points(X, R,
         pch = 16,
         col = cols[5])
  legend(
    'topleft', cex=0.8,
    legend = c('residuals', '95% Conf. Int.', '95% Pred. Int.'),
    pch = c(16, NA, NA),
    col = c(cols[5], cols_tr2[2:3]),
    lty = c(1, 1, 1),
    lwd = c(1, 10, 10),
    bty = 'o', box.col = NA, bg='white'
  )
  box()

  # # Coverage statistics
  # nSeg = length(cp$m)
  # ylim = c(min(cp$q95),1.0)
  # plot(1:nSeg,cp$m, ylim=ylim, col=4, pch=19,
  #      xlab = 'X-segment #',
  #      ylab = 'Probability',
  #      main = 'Coverage probability of 95% PUI-CI')
  # grid()
  # abline(h=0.95,col=2,lty=2)
  # segments(1:nSeg,cp$q95[,1],1:nSeg,cp$q95[,2],col=4)
}




