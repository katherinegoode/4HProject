newdata<-dicretereplicatesasattributes[,-(1:3)]
newdata
n_c <- function(data){
  n_c0 <- length(which(data[,1] ==0))
  n_c1 <- length(which(data[,1] ==1))
  return(c(n_c0,n_c1))
}

#returns a vector of P(c0), P(c1)
P_c <- function(s, t_c0,data){
  n_c0 <- n_c(data)[1]
  n_c1 <- n_c(data)[2]
  N = n_c0+n_c1
  P_c0 <- (n_c0 + s*t_c0)/(N + s)
  P_c1 <- (n_c1 + s*(1 - t_c0))/(N + s)
  return(c(P_c0,P_c1))
}

#this function returns a table n(ai, cj) where the 0-1 at the top is cj of if they have the disease or not
# ai is the attribute where i should start from 2
n_Ai <- function(i, data){
  tab <- table(data[,i], data[,1])
  return(tab)
}

P_ai <- function(i, attrib, data, c, s, tc, taic){
  n_c <- n_c(data)[c+1]
  n_aic <- n_Ai(i,data)[attrib,c+1]
  prob <- (n_aic +s*taic)/(n_c + s*tc)
  return(prob)
}
#produces a vector of n(ai, c) for the case given, where the case should just be a list of attributes ie 
#have removed the class.
n_aiccase <- function(data, case, c){
  alphai <- c(n_Ai(2,data)[as.matrix(case[1]+1),c+1])
  for (i in 2:length(case)){
    a <- n_Ai(i+1,data)[as.matrix(case[i]),c+1]
    alphai <- append(alphai,a)
  }
  return(alphai)
}


#All functions below are testing if the case c = 0 (no disease) credal dominates the case c = 1 (PD)
h1 <- function(data,case, x, s){
  alpha <- n_c(data)[1]
  beta <- n_c(data)[2]
  k <- dim(case)[2]
  alphai <- n_aiccase(data, case, 0)
  betai <- n_aiccase(data, case, 1)
  f <- (((beta + x)/(alpha +s -x))^(k-1)) * prod(alphai/(betai+x))
  return(f)
}
lnh1 <- function(data,case, x, s){
  alpha <- n_c(data)[1]
  beta <- n_c(data)[2]
  k <- dim(case)[2]
  alphai <- n_aiccase(data, case, 0)
  betai <- n_aiccase(data, case, 1)
  f <- (k-1)*(log(beta+x)-log(alpha +s -x)) + sum(log(alphai))- sum(log(betai+x))
  return(f)
}
dxlnh1 <- function(data,case, x, s){
  alpha <- n_c(data)[1]
  beta <- n_c(data)[2]
  k <- dim(case)[2]
  betai <- n_aiccase(data, case, 1)
  if(prod(betai) == 0){
    if(x ==  0){
      return(-10*100)
    }
    else{
      f <- (k-1)/(beta+x) + (k-1)/(alpha +s -x) - sum(1/(betai+x))
      return(f)
    }
  }
  else{
    f <- (k-1)/(beta+x) + (k-1)/(alpha +s -x) - sum(1/(betai+x))
    return(f)
  }
}
infh1 <- function(data, case, s){
  if(prod(n_aiccase(data, case, 0)) == 0){
    return(0)
  }
  if(dxlnh1(data, case, 0, s)>= 0){
    return(h(data, case, 0, s))
  }
  if(dxlnh1(data, case, s, s)<= 0){
    return(h(data, case, s, s))
  }
  else{
    lnhcase < function(x){
      return(lnh(data, case, x, s))
    }
    st <- newtonRaphson(lnhcase, 0)
    return(h1(data, case, st, s))
  }
}

#All functions below are testing if the case c = 1 (PD) credal dominates the case c = 0 (no disease)
h2 <- function(data,case, x, s){
  alpha <- n_c(data)[2]
  beta <- n_c(data)[1]
  k <- dim(case)[2]
  alphai <- n_aiccase(data, case, 1)
  betai <- n_aiccase(data, case, 0)
  f <- (((beta + x)/(alpha +s -x))^(k-1)) * prod(alphai/(betai+x))
  return(f)
}
lnh2 <- function(data,case, x, s){
  alpha <- n_c(data)[2]
  beta <- n_c(data)[1]
  k <- dim(case)[2]
  alphai <- n_aiccase(data, case, 1)
  betai <- n_aiccase(data, case, 0)
  f <- (k-1)*(log(beta+x)-log(alpha +s -x)) + sum(log(alphai))- sum(log(betai+x))
  return(f)
}
dxlnh2 <- function(data,case, x, s){
  alpha <- n_c(data)[2]
  beta <- n_c(data)[1]
  k <- dim(case)[2]
  betai <- n_aiccase(data, case, 0)
  if(prod(betai) == 0){
    if(x ==  0){
      return(-10*100)
    }
    else{
      f <- (k-1)/(beta+x) + (k-1)/(alpha +s -x) - sum(1/(betai+x))
      return(f)
    }
  }
  else{
    f <- (k-1)/(beta+x) + (k-1)/(alpha +s -x) - sum(1/(betai+x))
    return(f)
  }
}
infh2 <- function(data, case, s){
  if(prod(n_aiccase(data, case, 1)) == 0){
    return(0)
  }
  if(dxlnh2(data, case, 0, s)>= 0){
    return(h2(data, case, 0, s))
  }
  if(dxlnh2(data, case, s, s)<= 0){
    return(h2(data, case, s, s))
  }
  else{
    lnhcase < function(x){
      return(lnh2(data, case, x, s))
    }
    st <- newtonRaphson(lnhcase, 0)
    return(h2(data, case, st, s))
  }
}


#Now test if c0 dominates c1, if c1 dominates c0, or neither dominate either, and return a classification:

classificcredal <- function(data, case, s){
  if(infh1(data, case, s) > 1){
    return(c(0))
  }
  if(infh2(data,case, s) > 1){
    return(c(1))
  }
  else{
    return(c('indeterminate'))
  }
}


predccredal<- function(data, obs, s){
  p <- c()
  for (i in 1:dim(obs)[1]){
    c <- classificcredal(data, obs[i,],s)
    p <- append(p,c)
  }
  return(p)
}

#Uniform prior for NBC
g <- function(data, case){
  k <- length(case)
  numerator0 <- n_aiccase(data, case, 0)
  numerator0[1] <- numerator0[1] + 0.25
  numerator0 <- numerator0 + 0.25
  numerator1 <- n_aiccase(data, case, 1)
  numerator1[1] <- numerator1[1] + 0.25
  numerator1 <- numerator1 + 0.25
  g0 <- ((n_c(data)[1] + 1)^(1-k))*prod(numerator0)
  g1 <- ((n_c(data)[2] + 1)^(1-k))*prod(numerator1)
  return(g0 - g1)
}

classificbayes <- function(data, case){
  if(g(data, case) > 0){
    return(c(0))
  }
  else{
    return(c(1))
  }
}

predcbayes<- function(data, obs){
  p <- c()
  for (i in 1:dim(obs)[1]){
    c <- classificbayes(data, obs[i,])
    p <- append(p,c)
  }
  return(p)
}

tencross <- function(i, data, test){
  n = dim(data)[1]
  len <- n/10
  start <- (i-1)*len+1
  fin <- i*len
  if(test == 'training'){
    return(data[-(start:fin),])
  }
  if(test == 'test'){
    return(data[(start:fin),])
  }
}


crossvalidatecredal <- function(data, s){
  rdm <- sample(80,80, replace = FALSE)
  randomorder <- c()
  for (i in (1:80)){
    randomorder <- append(randomorder, rdm[i])
  }
  rdmdata <- data[randomorder,]
  classescredal <- c()
  classesbayes <- c()
  for (i in (1:10)){
    test <- tencross(i, rdmdata, 'test')[,-1]
    training <- tencross(i, rdmdata, 'training')
    classescredal <- append(classescredal, predccredal(training, test, s))
    classesbayes <- append(classesbayes, predcbayes(training, test))
  }
  comparisoncredal <- cbind(rdmdata[,1],classescredal)
  comparisonbayes <- cbind(rdmdata[,1],classesbayes)
  checkcredal <- c()
  for (i in(1:dim(rdmdata)[1])){
    if (comparisoncredal[i,1] == comparisoncredal[i,2]){
      checkcredal <- append(checkcredal, as.numeric(1))
    }
    else{
      if (comparisoncredal[i,2] == 'indeterminate'){
        checkcredal <- append(checkcredal, as.numeric(-100))
      }
      else{
        checkcredal <- append(checkcredal, as.numeric(0))
      }
    }
  }
  checkbayes <- c()
  for (i in(1:dim(rdmdata)[1])){
    if (comparisonbayes[i,1] == comparisonbayes[i,2]){
      checkbayes <- append(checkbayes, as.numeric(1))
    }
    else{
      checkbayes <- append(checkbayes, as.numeric(0))
    }
  }
  resultcredal <- as.data.frame(cbind(randomorder, comparisoncredal, checkcredal, comparisonbayes[,2], checkbayes))
  ordres<- resultcredal[order(randomorder),] 
  ordresult <- as.data.frame(cbind(newdata[,1], ordres))
  names(ordresult)[1] <- 'PatientId'
  names(ordresult)[2] <- 'RecordingNumber'
  names(ordresult)[3] <- 'ActualClass'
  names(ordresult)[4] <- 'NCCPredictedClass'
  names(ordresult)[5] <- 'NCCAccurate?'
  names(ordresult)[6] <- 'NBCPredictedClass'
  names(ordresult)[7] <- 'NBCAccurate?'
  return(ordresult)
}


repasapp1 <- crossvalidatecredal(newdata,1)
com <- function(ordresult){
  indet <- ordresult[which(ordresult[,5] == -100), ]
  corr <- indet[which(indet[, 7] == 1),]
  numindet <- dim(indet)[1]
  numcorr <- dim(corr)[1]
  return(numcorr/numindet)
}

nccvcnbc <- c()
for(i in 1:6){
  nccvcnbc <- append(nccvcnbc, com(crossvalidatecredal(newdata, 1)))
}
length(nccvcnbc)
mean(nccvcnbc)

dim(repasapp1[which(repasapp1[,4] == 'indeterminate'), ])



accuracy1 <- function(predtable){
  correct <- sum(predtable[,5] == 1)
  indet <- sum(predtable[,5] == -100)
  N <- dim(predtable)[1]
  da <- (1/N)*(correct + 0.5*indet) 
  return(c(correct/(N-indet), 1-indet/N, da))
}


falsenegpos <- function(predtable){
  incorr <- predtable[which(predtable[,5] == 0),]
  fpos <- incorr[which(incorr[,3]==0),]
  fneg <- incorr[which(incorr[,3]==1),]
  Num <- dim(incorr)[1]
  p<- dim(fpos)[1]/Num
  n<- dim(fneg)[1]/Num
  return(c(p,n))
  
}
averageacc1 <- function(data, s, r){
  sing <- c()
  determ <- c()
  dacc <- c()
  falsepos <- c()
  falseneg <- c()
  for(i in 1:r){
    predtables <- crossvalidate(data, s)
    resu <- accuracy1(predtables)
    sing <- append(sing, resu[1])
    determ <- append(determ, resu[2])
    dacc <- append(dacc, resu[3])
    
    fal <- falsenegpos(predtables)
    falsepos <- append(falsepos, fal[1])
    falseneg <- append(falseneg, fal[2])
  }
  return(c(mean(sing),mean(determ), mean(dacc), mean(falsepos), mean(falseneg)))
}
averageacc1(newdata, 1, 20)
