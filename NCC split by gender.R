mydata <-naremoved[,-(1:2)]
menmydata <- mydata[which(mydata[,2]==0),-2]
womanmydata <- mydata[which(mydata[,2]==1),-2]
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
  alphai <- c()
  for (i in 1:length(case)){
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

classific <- function(data, case, s){
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

predc<- function(data, obs, s){
  p <- c()
  for (i in 1:dim(obs)[1]){
    c <- classific(data, obs[i,],s)
    p <- append(p,c)
  }
  return(p)
}

tencross <- function(i, data, test){
  n = dim(data)[1]
  len <- n/8
  start <- (i-1)*len+1
  fin <- i*len
  if(test == 'training'){
    return(data[-(start:fin),])
  }
  if(test == 'test'){
    return(data[(start:fin),])
  }
}

crossvalidate <- function(data, s){
  rdm <- sample(32,32, replace = FALSE)*3
  randomorder <- c()
  for (i in (1:32)){
    randomorder <- append(randomorder, c(rdm[i]-2,rdm[i]-1,rdm[i]))
  }
  rdmdata <- data[randomorder,]
  classes <- c()
  for (i in (1:8)){
    test <- tencross(i, rdmdata, 'test')[,-1]
    training <- tencross(i, rdmdata, 'training')
    classes <- append(classes, predc(training, test, s))
  }
  comparison <- cbind(rdmdata[,1],classes)
  check <- c()
  for (i in(1:dim(rdmdata)[1])){
    if (comparison[i,1] == comparison[i,2]){
      check <- append(check, as.numeric(1))
    }
    else{
      if (comparison[i,2] == 'indeterminate'){
        check <- append(check, as.numeric(-100))
      }
      else{
        check <- append(check, as.numeric(0))
      }
    }
  }
  result <- as.data.frame(cbind(randomorder, comparison, check))
  ordres<- result[order(randomorder),] 
  ordresult <- as.data.frame(cbind(seq(1,32,1), ordres))
  names(ordresult)[1] <- 'PatientId'
  names(ordresult)[2] <- 'RecordingNumber'
  names(ordresult)[3] <- 'ActualClass'
  names(ordresult)[4] <- 'PredictedClass'
  names(ordresult)[5] <- 'Accurate?'
  return(ordresult)
}

patientprediction <- function(prediction){
  cfic <- c()
  for (i in(1:32)){
    cfic <- append(cfic, sum(as.numeric(as.character(prediction$`Accurate?`[(3*i-2):(3*i)]))))
  }
  for (i in (1:32)){
    if(cfic[i] == 3){
      cfic[i] = 'correct'
    }
    if(cfic[i] == 2){
      cfic[i] = 'correct'
    }
    if(cfic[i] == 1){
      cfic[i] = 'wrong'
    }
    if(cfic[i] == 0){
      cfic[i] = 'wrong'
    }
    if(cfic[i] == -98){
      cfic[i] = 'correct'
    }
    if(cfic[i] == -99){
      cfic[i] = 'indeterminate'
    }
    if(cfic[i] == -100){
      cfic[i] = 'wrong'
    }
    if(cfic[i] == -300){
      cfic[i] = 'indeterminate'
    }
  }
  actclass <- c()
  for (i in(1:32)){
    actclass <- append(actclass, prediction$ActualClass[3*i])
  }
  predictions <- as.data.frame(cbind(actclass, cfic))
  names(predictions)[1] <- 'class'
  names(predictions)[2] <- 'PredictionAccuracy'
  return(predictions)
}

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

accuracy2 <- function(predtable){
  correct <- sum(predtable[,2] == 'correct')
  indet <- sum(predtable[,2] == 'indeterminate')
  N <- dim(predtable)[1]
  da <- (1/N)*(correct + 0.5*indet)
  return(c(correct/(N-indet), 1-indet/N, da))
}

falsenegpos2 <- function(predtable){
  incorr <- predtable[which(predtable[,2] == 'wrong'),]
  fpos <- incorr[which(incorr[,1]==0),]
  fneg <- incorr[which(incorr[,1]==1),]
  Num <- dim(incorr)[1]
  p<- dim(fpos)[1]/Num
  n<- dim(fneg)[1]/Num
  return(c(p,n))
}

averageacc2 <- function(data, s, r){
  sing <- c()
  determ <- c()
  dacc <- c()
  falsepos <- c()
  falseneg <- c()
  for(i in 1:r){
    predtables <- patientprediction(crossvalidate(data, s))
    resu <- accuracy2(predtables)
    sing <- append(sing, resu[1])
    determ <- append(determ, resu[2])
    dacc <- append(dacc, resu[3])
    
    fal <- falsenegpos2(predtables)
    falsepos <- append(falsepos, fal[1])
    falseneg <- append(falseneg, fal[2])
  }
  return(c(mean(sing),mean(determ), mean(dacc), mean(falsepos), mean(falseneg)))
}
averageacc1(mydata, 1, 20)
averageacc2(mydata, 1, 20)

a<- c(c(1,2), c(3,4))
a

b <- c(averageacc1(menmydata, 1, 20), averageacc2(menmydata, 1, 20))

b
c <- c(averageacc1(womanmydata, 1, 20), averageacc2(womanmydata, 1, 20))
c
