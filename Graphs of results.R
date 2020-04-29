# ncc with s=1

singaccncc <- c(81.24, 85.93, 82.49)
determncc <- c(97.62, 97.69, 70.31)
discaccncc <- c(80.5, 85.09, 72.84)

counts <- cbind(singaccncc, determncc, discaccncc)
barplot(counts, main="NCC with s=1",
        ylab="%",ylim = c(0,100), 
        xlab = "Diagnosic Measure", 
        names.arg = c("Single Acc.", "Determ.", "Disc. Acc"), col=c("darkblue","red", "black"), 
        beside=TRUE, 
        legend = c(1,2,3),
        args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))




#compar
nbcaccur <- c(80.6, 85.25, 84.75)
nccdiscacc <- c(80.5, 85.09, 72.84)
indetacc <- c(46.97, 53.04, 89.16)
counts <- cbind(nbcaccur, nccdiscacc, indetacc)
barplot(counts, main="NBC v.s. NCC",
        ylab="%",ylim = c(0,100), 
        xlab = "Diagnosic Measure", 
        names.arg = c("NBC Acc.", "NCC D.Acc", "Indet. Acc"), col=c("darkblue","red", "black"), 
        beside=TRUE, 
        legend = c(1,2,3),
        args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))



overallres <- c(85.25, 85.09, 73.60, 75.2)
mensres <- c(81.77, 82.67, 68.37, 70.6)
womensres <- c(87.81, 82.19, 80.68, 87.6)
plot(c(1,2,3,4),overallres, xaxt = "n", main = "Comparing Classifiers", ylab="%", ylim = c(68.2, 87.9), 
     xlab = "Classifier",
     col= "black", pch = 16)
axis(1, at = seq(1, 4, by = 1), labels =c("NBC", "NCC", "MLE", "GLM"), las=1)
points(seq(4), mensres, col = "red", pch = 16)
points(seq(4), womensres, col = "green", pch = 16)

legend(x=1,y=78,c("All Data", "Men", "Women"),cex=.8,col=c("black","red", "green"),pch=c(16,16))

min(overallres, mensres, womensres)

plot(seq(40), ParkinsonData$Jitter_rel[3*c(1:40)], xlim = c(0,80))
points(seq(40), ParkinsonData$Jitter_rel[3*c(1:40)-1])
points(seq(40), ParkinsonData$Jitter_rel[3*c(1:40)-2])
points(40 + seq(40), ParkinsonData$Jitter_rel[3*c(41:80)], col = "red")
points(40 + seq(40), ParkinsonData$Jitter_rel[3*c(41:80)-1], col = "red")
points(40 + seq(40), ParkinsonData$Jitter_rel[3*c(41:80)-2], col = "red")




cont <- ParkinsonData$Jitter_rel[c(1:120)]
park <- ParkinsonData$Jitter_rel[c(121:240)]
x <- list("CONT"=cont, "PARK"=park)
stripchart(x, xlim = c(0,2.6),
           main="Relative Jitter Comparison",
           xlab="Relative Jitter",
           ylab = "Patient",
           method="jitter",
           col=c("orange","red"),
           pch=16
)
cont <- ParkinsonData$Shim_loc[c(1:120)]
park <- ParkinsonData$Shim_loc[c(121:240)]
x <- list("CONT"=cont, "PARK"=park)
stripchart(x,
           main="Local Shimmer Comparison",
           xlab="Local Shimmer",
           ylab = "Patient",
           method="jitter",
           col=c("orange","red"),
           pch=16
)

cont <- ParkinsonData$HNR38[c(1:120)]
park <- ParkinsonData$HNR38[c(121:240)]
x <- list("CONT"=cont, "PARK"=park)
stripchart(x,
           main="Harmonic-to-Noise Ratio Comparison",
           xlab="Harmonic-to-Noise Ratio",
           ylab = "Patient",
           method="jitter",
           col=c("orange","red"),
           pch=16
)
cont <- ParkinsonData$MFCC0[c(1:120)]
park <- ParkinsonData$MFCC0[c(121:240)]
x <- list("CONT"=cont, "PARK"=park)
stripchart(x,
           main="MFCC Comparison",
           xlab="MFCC",
           ylab = "Patient",
           method="jitter",
           col=c("orange","red"),
           pch=16
)
cont <- ParkinsonData$Delta12[c(1:120)]
park <- ParkinsonData$Delta12[c(121:240)]
x <- list("CONT"=cont, "PARK"=park)
stripchart(x,
           main="Delta Comparison",
           xlab="Delta",
           ylab = "Patient",
           method="jitter",
           col=c("orange","red"),
           pch=16
)




plot(c(1,2,3), c(80.38, 84.94, 84.44), xaxt = "n", ylim = c(80.3, 85.3), pch = 16, xlab = "Method", 
     ylab = "Accuracy", main = "Jeffrey's v.s Unform prior")
axis(1, at = seq(1, 3, by = 1), las=1)
points(c(80.60, 85.25, 84.75), col = "red", pch = 16)

legend(x=1,y=85,c("Jeffrey's", "Uniform"),cex=.8,col=c("black","red"),pch=c(16,16))




?legend()
points(c(80.52, 84.80, 72), col = "blue", pch = 16)


