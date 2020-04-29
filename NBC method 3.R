newdata<-dicretereplicatesasattributes[,-(1:3)]
mennewdata <- newdata[which(newdata[,2]==0),-2]
womannewdata <- newdata[which(newdata[,2]==1),-2]


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
  numerator0[1] <- numerator0[1] + 0.25
  numerator0 <- numerator0 + 0.25
  numerator1 <- n_aiccase(data, case, 1)
  numerator1[1] <- numerator1[1] + 0.25
  numerator1 <- numerator1 + 0.25
  g0 <- ((n_c(data)[1] + 1)^(1-k))*prod(numerator0)
  g1 <- ((n_c(data)[2] + 1)^(1-k))*prod(numerator1)
  return(g0 - g1)
}


classific <- function(data, case){
  if(g(data, case) > 0){
    return(c(0))
  }
  else{
    return(c(1))
  }
}

predc<- function(data, obs){
  p <- c()
  for (i in 1:dim(obs)[1]){
    c <- classific(data, obs[i,])
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

crossvalidate <- function(data){
  rdm <- sample(32,32, replace = FALSE)
  randomorder <- c()
  for (i in (1:32)){
    randomorder <- append(randomorder, rdm[i])
  }
  rdmdata <- data[randomorder,]
  classes <- c()
  for (i in (1:8)){
    test <- tencross(i, rdmdata, 'test')[,-1]
    training <- tencross(i, rdmdata, 'training')
    classes <- append(classes, predc(training, test))
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
  ordresult <- as.data.frame(cbind(seq(1,32,1), ordres))
  names(ordresult)[1] <- 'PatientId'
  names(ordresult)[2] <- 'RecordingNumber'
  names(ordresult)[3] <- 'ActualClass'
  names(ordresult)[4] <- 'PredictedClass'
  names(ordresult)[5] <- 'Accurate?'
  return(ordresult)
}


accuracy <- function(predtable){
  correct <- sum(predtable[,5] == 1)
  indet <- sum(predtable[,5] == -100)
  N <- dim(predtable)[1]
  return((correct/(N-indet))*100)
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
averageacc1 <- function(data){
  p <- c()
  fp <- c()
  fn <- c()
  for(i in 1:20){
    predtables <-crossvalidate(data)
    fpresu <- falsenegpos(predtables)
    p <- append(p, accuracy(predtables))
    fp <- append(fp, fpresu[1])
    fn <- append(fn, fpresu[2])
  }
  return(c(mean(p), mean(fp), mean(fn)))
}

averageacc1(womannewdata)













