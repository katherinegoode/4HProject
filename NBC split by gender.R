mydata <-naremoved[,-(1:2)]
menmydata <- mydata[which(mydata[,2]==0),-2]
womanmydata <- mydata[which(mydata[,2]==1),-2]

n_c <- function(data){
  n_c0 <- length(which(data[,1] ==0))
  n_c1 <- length(which(data[,1] ==1))
  return(c(n_c0,n_c1))
}

#this function returns a table n(ai, cj) where the 0-1 at the top is cj of if they have the disease or not
# ai is the attribute where i should start from 2
n_Ai <- function(i, data){
  tab <- table(data[,i], data[,1])
  return(tab)
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

#uniform prior
g <- function(data, case){
  k <- length(case)
  numerator0 <- n_aiccase(data, case, 0)
  numerator0 <- numerator0 + 0.25
  numerator1 <- n_aiccase(data, case, 1)
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

dim(womanmydata)
96/3
crossvalidatebayes <- function(data){
  rdm <- sample(48,48, replace = FALSE)*3
  randomorder <- c()
  for (i in (1:48)){
    randomorder <- append(randomorder, c(rdm[i]-2,rdm[i]-1,rdm[i]))
  }
  rdmdata <- data[randomorder,]
  classes <- c()
  for (i in (1:8)){
    test <- tencross(i, rdmdata, 'test')[,-1]
    training <- tencross(i, rdmdata, 'training')
    classes <- append(classes, predcbayes(training, test))
  }
  comparison <- cbind(rdmdata[,1],classes)
  check <- c()
  for (i in(1:dim(rdmdata)[1])){
    if (comparison[i,1] == comparison[i,2]){
      check <- append(check, as.numeric(1))
    }
    else{
      check <- append(check, as.numeric(0))
    }
  }
  result <- as.data.frame(cbind(randomorder, comparison, check))
  ordres<- result[order(randomorder),] 
  ordresult <- as.data.frame(cbind(seq(1,48,1), ordres))
  names(ordresult)[1] <- 'PatientId'
  names(ordresult)[2] <- 'RecordingNumber'
  names(ordresult)[3] <- 'ActualClass'
  names(ordresult)[4] <- 'BayesPredictedClass'
  names(ordresult)[5] <- 'Accurate?'
  return(ordresult)
}
crossvalidatebayes(womanmydata)
patientprediction <- function(prediction){
  cfic <- c()
  for (i in(1:48)){
    cfic <- append(cfic, sum(as.numeric(as.character(prediction$`Accurate?`[(3*i-2):(3*i)]))))
  }
  for (i in (1:48)){
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
  }
  actclass <- c()
  for (i in(1:48)){
    actclass <- append(actclass, prediction$ActualClass[3*i])
  }
  predictions <- as.data.frame(cbind(actclass, cfic))
  names(predictions)[1] <- 'class'
  names(predictions)[2] <- 'PredictionAccuracy'
  return(predictions)
}


accuracy1 <- function(predtable){
  correct <- sum(predtable[,5] == 1)
  N <- dim(predtable)[1]
  return((correct/N)*100)
}

accuracy2 <- function(predtable){
  correct <- sum(predtable[,2] == 'correct')
  N <- dim(predtable)[1]
  return((correct/N)*100)
}
falsenegpos1 <- function(predtable){
  incorr <- predtable[which(predtable[,5] == 0),]
  fpos <- incorr[which(incorr[,3]==0),]
  fneg <- incorr[which(incorr[,3]==1),]
  Num <- dim(incorr)[1]
  p<- dim(fpos)[1]/Num
  n<- dim(fneg)[1]/Num
  return(c(p,n))
  
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

averageacc1 <- function(data){
  p <- c()
  fp <- c()
  fn <- c()
  for(i in 1:20){
    predtables <-crossvalidatebayes(data)
    fpresu <- falsenegpos1(predtables)
    p <- append(p, accuracy1(predtables))
    fp <- append(fp, fpresu[1])
    fn <- append(fn, fpresu[2])
  }
  return(c(mean(p), mean(fp), mean(fn)))
}


averageacc2 <- function(data){
  p <- c()
  fp <- c()
  fn <- c()
  for(i in 1:20){
    predtables <-patientprediction(crossvalidatebayes(data)) 
    fpresu <- falsenegpos2(predtables)
    p <- append(p, accuracy2(predtables))
    fp <- append(fp, fpresu[1])
    fn <- append(fn, fpresu[2])
  }
  return(c(mean(p), mean(fp), mean(fn)))
}


resu2 <- c(averageacc1(menmydata),averageacc2(menmydata))
resu2



