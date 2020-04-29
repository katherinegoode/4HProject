pddata <- ParkinsonData[,-(1:2)]

menpddata <- pddata[which(pddata[,2]==0),-2]
womanpddata <- pddata[which(pddata[,2]==1),-2]


mu <- function(data, k){
  neg <- c()
  pos <- c()
  for(i in 1:dim(data)[1]){
    if(data[i,1] == 0){
      neg <- append(neg, i)
    }
    else{
      pos <- append(pos, i)
    }
  }
  cont <- data[neg,]
  park <- data[pos,]
  mucont <- sum(cont[,k])/length(cont[,k])
  mupark <- sum(park[,k])/length(park[,k])
  varcont<- sum((cont[,k] - mucont)^2)/length(cont[,k])
  varpark<- sum((park[,k] - mupark)^2)/length(park[,k])
  return(c(mucont, varcont, mupark, varpark))
}

lambdafn <- function(data, case){
  ex <- c()
  pro <- c()
  for (i in 2:dim(data)[2]){
    val <- mu(data, i)
    pro <- append(pro, (val[4]/val[2])^2)
    bigfrac <- ((case[i-1] - val[3])^2)/(val[4]^0.5) - ((case[i-1] - val[1])^2)/(val[2]^0.5)
    ex <- append(ex, bigfrac)
  }
  lam <- prod(pro)*exp(0.5*sum(ex))
  return(lam)
}

predc<- function(data, obs){
  p <- c()
  for (i in 1:dim(obs)[1]){
    c <- lambdafn(data, obs[i,])
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



lamdastar <- function(pred, lamdast){
  predic1 <- c()
  for (i in 1:96){
    if (pred[i,2] > lamdast){
      predic1 <- append(predic1, 0)
    }
    if (pred[i,2] < 1/lamdast){
      predic1 <- append(predic1, 1)
    }
    if(1/lamdast < pred[i,2]){
      if (pred[i,2] <lamdast){
        predic1 <- append(predic1, -100)
      }
    }
  }
  return(predic1)
}


crossvalidate <- function(data, lamdast){
  rdm <- sample(32, 32, replace = FALSE)*3
  randomorder <- c()
  for (i in (1:32)){
    randomorder <- append(randomorder, c(rdm[i]-2,rdm[i]-1,rdm[i]))
  }
  rdmdata <- data[randomorder,]
  lambd <- c()
  for (i in (1:8)){
    test <- tencross(i, rdmdata, 'test')[,-1]
    training <- tencross(i, rdmdata, 'training')
    lambd <- append(lambd, predc(training, test))
  }
  comparison <- cbind(rdmdata[,1],lambd)
  result <- as.data.frame(cbind(randomorder, comparison))
  ordres<- result[order(randomorder),] 
  ordres <- ordres[, -1]
  names(ordres)[1] <- 'ActualClass'
  names(ordres)[2] <- 'Lamba'
  accurc <- cbind(ordres, lamdastar(ordres, lamdast))
  names(accurc)[3] <- 'Predicted Class'
  check <- c()
  for (i in(1:dim(accurc)[1])){
    if (accurc[i,1] == accurc[i,3]){
      check <- append(check, 'correct')
    }
    else{
      if (accurc[i,3] == -100){
        check <- append(check, 'indeterminate')
      }
      else{
        check <- append(check, 'wrong')
      }
    }
  }
  final <- cbind(accurc, check)
  names(final)[4] <- 'correct?'
  return(final)
}

accuracy <- function(predtable){
  correct <- sum(predtable[,4] == 'correct')
  indet <- sum(predtable[,4] == 'indeterminate')
  N <- dim(predtable)[1]
  da <- (1/N)*(correct + 0.5*indet) 
  return(c(correct/(N-indet), 1-indet/N, da))
}

falsenegpos <- function(predtable){
  incorr <- predtable[which(predtable[,4] == 'wrong'),]
  fpos <- incorr[which(incorr[,1]==0),]
  fneg <- incorr[which(incorr[,1]==1),]
  Num <- dim(incorr)[1]
  p<- dim(fpos)[1]/Num
  n<- dim(fneg)[1]/Num
  return(c(p,n))
  
}


averageacc <- function(data,lamdast, r){
  sing <- c()
  determ <- c()
  dacc <- c()
  falsepos <- c()
  falseneg <- c()
  for(i in 1:r){
    predtables <- crossvalidate(data, lamdast)
    resu <- accuracy(predtables)
    sing <- append(sing, resu[1])
    determ <- append(determ, resu[2])
    dacc <- append(dacc, resu[3])
    fal <- falsenegpos(predtables)
    falsepos <- append(falsepos, fal[1])
    falseneg <- append(falseneg, fal[2])
  }
  return(c(mean(sing),mean(determ), mean(dacc), mean(falsepos), mean(falseneg)))
}
averageacc(womanpddata, 1000, 20)


s3 <- c()
d3 <-c()
da3 <-c()
for(i in c(193:200)){
  resu <- accuracy(crossvalidate(pddata, i))
  s3 <- append(s3, resu[1])
  d3 <- append(d3, resu[2])
  da3 <- append(da3, resu[3])
}

predictions <- function(results, lamdast){
  accurc <- cbind(results, lamdastar(results, lamdast))
  names(accurc)[3] <- 'Predicted Class'
  check <- c()
  for (i in(1:dim(accurc)[1])){
    if (accurc[i,1] == accurc[i,3]){
      check <- append(check, 'correct')
    }
    else{
      if (accurc[i,3] == -100){
        check <- append(check, 'indeterminate')
      }
      else{
        check <- append(check, 'wrong')
      }
    }
  }
  final <- cbind(accurc, check)
  names(final)[4] <- 'correct?'
  return(final)
}


for (j in c(1:5)){
  testing <- crossvalidate(pddata, 2)
  s <- c()
  d <- c()
  da<- c()
  for(i in c(1:2000)){
    resu <- accuracy(predictions(testing[,(1:2)], i))
    s <- append(s, resu[1])
    d <- append(d, resu[2])
    da <- append(da, resu[3])
  }
  singavv <- singavv + s
  determavv <- determavv + d
  daccavv <- daccavv + da
}
par(mfrow = c(3,3))

determavv

singavv20[c(500, 1000, 1500, 2000)]
determavv20[c(500, 1000, 1500, 2000)]*100
daccavv20[c(500, 1000, 1500, 2000)]*100

max(singavv5)
ylim =c(0.73, 0.7366667)
plot(c(1:2000),singavv5, type = 'l', xlab = 'lamda*', ylab='Single Acc.', main = 'Single Accuracy: r = 5', ylim =c(0.7325, 0.748054))
plot(c(1:2000),determavv5, type = 'l', xlab = 'lamda*', ylab='Determinacy', main = 'Determinacy: r = 5', ylim =c(0.95, 1))
plot(c(1:2000),daccavv5, type = 'l', xlab = 'lamda*', ylab='d-acc.', main = 'Discounted Accuracy: r = 5', ylim =c(0.73, 0.7366667))

plot(c(1:2000),singavv10, type = 'l', xlab = 'lamda*', ylab='Single Acc.', main = 'Single Accuracy: r=10', ylim =c(0.7325, 0.748054))
plot(c(1:2000),determavv10, type = 'l', xlab = 'lamda*', ylab='Determinacy', main = 'Determinacy: r = 10', ylim =c(0.95, 1))
plot(c(1:2000),daccavv10, type = 'l', xlab = 'lamda*', ylab='d-acc.', main = 'Discounted Accuracy: r = 10', ylim =c(0.73, 0.7366667))



plot(c(1:2000),singavv20, type = 'l', xlab = 'lamda*', ylab='Single Acc.', main = 'Single Accuracy: r = 20', ylim =c(0.7325, 0.748054))
plot(c(1:2000),determavv20, type = 'l', xlab = 'lamda*', ylab='Determinacy', main = 'Determinacy: r =20', ylim =c(0.95, 1))
plot(c(1:2000),daccavv20, type = 'l', xlab = 'lamda*', ylab='d-acc.', main = 'Discounted Accuracy: r= 20', ylim =c(0.73, 0.7366667))

max(daccavv5, daccavv10, daccavv20)
