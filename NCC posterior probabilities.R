mydata <-naremoved[,-(1:2)]

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


h1 <- function(data,case, x, s){
  alpha <- n_c(data)[1]
  beta <- n_c(data)[2]
  k <- dim(case)[2]
  alphai <- n_aiccase(data, case, 0)
  betai <- n_aiccase(data, case, 1)
  if(prod(alpha*beta*betai*alphai) == 0){
    return(0)
  }
  else{
    f <- (((beta + x)/(alpha +s -x))^(k-1)) * prod(alphai/(betai+x))
    return(f)
  }
}

h2 <- function(data,case, x, s){
  alpha <- n_c(data)[2]
  beta <- n_c(data)[1]
  k <- dim(case)[2]
  alphai <- n_aiccase(data, case, 1)
  betai <- n_aiccase(data, case, 0)
  if(prod(alpha*beta*betai*alphai) == 0){
    return(0)
  }
  else{
    f <- (((beta + x)/(alpha +s -x))^(k-1)) * prod(alphai/(betai+x))
    return(f)
  }
}

#Now test if c0 dominates c1, if c1 dominates c0, or neither dominate either, and return a classification:
classific <- function(data, case, s){
  h1case<- function(x){
    return(h1(data,case, x, s))
  }
  h2case<- function(x){
    return(h2(data,case, x, s))
  }
  if(optimise(h1case, interval=c(0, s))$objective > 1){
    return(c(0))
  }
  if(optimise(h2case, interval=c(0, s))$objective > 1){
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

crossvalidate <- function(data, s){
  rdm <- sample(80,80, replace = FALSE)*3
  randomorder <- c()
  for (i in (1:80)){
    randomorder <- append(randomorder, c(rdm[i]-2,rdm[i]-1,rdm[i]))
  }
  rdmdata <- data[randomorder,]
  classes <- c()
  for (i in (1:10)){
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
  ordresult <- as.data.frame(cbind(naremoved[,1], ordres))
  names(ordresult)[1] <- 'PatientId'
  names(ordresult)[2] <- 'RecordingNumber'
  names(ordresult)[3] <- 'ActualClass'
  names(ordresult)[4] <- 'PredictedClass'
  names(ordresult)[5] <- 'Accurate?'
  return(ordresult)
}

patientprediction <- function(prediction){
  cfic <- c()
  for (i in(1:80)){
    cfic <- append(cfic, sum(as.numeric(as.character(prediction$`Accurate?`[(3*i-2):(3*i)]))))
  }
  for (i in (1:80)){
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
  predictions <- as.data.frame(cbind(levels(naremoved$ID), cfic))
  names(predictions)[1] <- 'PatientId'
  names(predictions)[2] <- 'PredictionAccuracy'
  return(predictions)
}

accuracy1 <- function(predtable){
  correct <- sum(predtable[,5] == 1)
  indet <- sum(predtable[,5] == -100)
  N <- dim(predtable)[1]
  return((correct/(N-indet))*100)
}

accuracy2 <- function(predtable){
  correct <- sum(predtable[,2] == 'correct')
  N <- dim(predtable)[1]
  return((correct/N)*100)
}

averageacc1 <- function(data, s){
  p <- c()
  for(i in 1:20){
    p <- append(p, accuracy1(crossvalidate(data, s)))
  }
  return(sum(p)/length(p))
}

averageacc2 <- function(data, s){
  p <- c()
  for(i in 1:20){
    p <- append(p, accuracy2(patientprediction(crossvalidate(data, s))))
  }
  return(sum(p)/length(p))
}

averageacc1(mydata, 1)
averageacc2(mydata, 1)
