naremoved
plot(continuousdata$Gender[-173], continuousdata$Jitter_RAP[-173])
continuousdata <- ReplicatedAcousticFeatures_ParkinsonDatabase_1_
identify(continuousdata$Gender, continuousdata$Jitter_RAP)
table(naremoved$Status, naremoved$Jitter_rel)
chisq.test(naremoved$Gender[(121:240)], naremoved$Jitter_rel[(121:240)])
