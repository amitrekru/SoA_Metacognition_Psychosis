library(stringr)
library(plyr)
library(readr)
library(tidyverse)
library(roperators)
library(raster)
library(foreach)

BetaClassifier <- function(LeaveOutNum, Iterations, fraction, scrambled, expFile, seed)
{
  if (missing(scrambled))
    scrambled = F
  
  if (missing(fraction))
    fraction = 1
  
  if (missing(expFile))
    fullExpFile <- read.csv("../Data/All_Participants_Trial_By_Trial_Data.csv", TRUE)
  else
    fullExpFile = expFile
  
  if (missing(seed))
    seed = 7
  
  set.seed(seed)
  
  print(sprintf("Classifier: Leave out %.d, %.d iterations, %g%% of each trial kind. Current time: %s", LeaveOutNum, Iterations, fraction*100, format(Sys.time(), "%X")))
  
  LeaveOutCount = LeaveOutNum
  IterationsCount = Iterations
  
  # Counters
  CorrectHc     = 0 # True negative
  CorrectScz    = 0 # True positive
  IncorrectHc   = 0 # False positive
  IncorrectScz  = 0 # False negative
  
  # Loading data
  fullExpFile$group<-ifelse(str_detect(fullExpFile$sub_name, 'Amit_HC'),"Healthy","Psychosis")
  fullExpFile$domain_name<-ifelse(fullExpFile$domain == 1, "Temporal","Spatial")
  fullExpFile$sub_name = str_c(fullExpFile$group, numextract(fullExpFile$sub_name), sep = " ")
  
  
  randHealthyId = sample(unique(fullExpFile[fullExpFile$group == "Healthy", "sub_name"]), 15);
  randPsychosisId = sample(unique(fullExpFile[fullExpFile$group == "Psychosis", "sub_name"]), 15);
  
  # Creating an empty dataframe of correct/incorrect counters
  sub_vector <- data.frame(matrix(ncol = 3, nrow = length(unique(fullExpFile$sub_name))))
  cols <- c("sub_name", "Correct", "Incorrect")
  colnames(sub_vector) <- cols
  sub_vector$sub_name = unique(fullExpFile$sub_name)
  sub_vector$Correct = 0
  sub_vector$Incorrect = 0
  
  # Creating an empty dataframe of correct/incorrect counters for each iteration
  template_iteration_sub_vector <- sub_vector
  
  iterationAccuracy = array(dim = IterationsCount)
  
  # Calculating mean SoA of every subject in both domains in all magnitudes
  IndividualSoA <- ddply(fullExpFile, .(sub_name,domain_name,mag), summarise,
                         SoA = mean(response)
  )
  
  sub_lm <- data.frame(matrix(ncol = 3, nrow = length(unique(fullExpFile$sub_name))))
  cols <- c("sub_name", "Temporal_lm_slope", "Spatial_lm_slope")
  colnames(sub_lm) <- cols
  sub_lm$sub_name = unique(fullExpFile$sub_name)
  
  for (subIndex in 1:nrow(sub_lm))
  {
    sub_lm[subIndex, "Temporal_lm_slope"] = lm(formula = SoA~mag, data = IndividualSoA[IndividualSoA$sub_name ==                    sub_lm[subIndex,"sub_name"] & IndividualSoA$domain_name == "Temporal",])$coefficients[2]
    sub_lm[subIndex, "Spatial_lm_slope"] = lm(formula = SoA~mag, data = IndividualSoA[IndividualSoA$sub_name ==                    sub_lm[subIndex,"sub_name"] & IndividualSoA$domain_name == "Spatial",])$coefficients[2]
  }
  
  
  # Running IterationsCount iterations
  for (IterationIndex in 1:IterationsCount)
  {
    iteration_sub_vector <- template_iteration_sub_vector
    
    expFile <- fullExpFile %>% group_by(sub_name, domain, mag) %>% sample_frac(fraction)
    
    # Changing to more convenient strings
    IndividualSoA$group = ifelse(str_detect(IndividualSoA$sub_name, 'Healthy'),"Healthy","Psychosis")
    
    if (scrambled)
    {
      IndividualSoA[which(IndividualSoA$sub_name %in% randHealthyId), "group"] = "Psychosis"
      IndividualSoA[which(IndividualSoA$sub_name %in% randPsychosisId), "group"] = "Healthy"
    }
    
    # Sampling LeaveOutCount participants from each group 
    HcIds = sample(x = unique(IndividualSoA[IndividualSoA$group == "Healthy", "sub_name"]), size = LeaveOutCount)
    SczIds = sample(x = unique(IndividualSoA[IndividualSoA$group == "Psychosis", "sub_name"]), size = LeaveOutCount)
    
    # Removing the sampled participants from the group dataset
    HcSubset = IndividualSoA[IndividualSoA$group == "Healthy" & !IndividualSoA$sub_name %in% HcIds,]
    SczSubset = IndividualSoA[IndividualSoA$group == "Psychosis" & !IndividualSoA$sub_name %in% SczIds,]
    
    # Computing group linear model of SoA for both groups for both domains separately
    HcTemporalModel = lm(formula = SoA~mag, data = HcSubset[HcSubset$domain_name == "Temporal",])
    HcSpatialModel = lm(formula = SoA~mag, data = HcSubset[HcSubset$domain_name == "Spatial",])
    SczTemporalModel = lm(formula = SoA~mag, data = SczSubset[SczSubset$domain_name == "Temporal",])
    SczSpatialModel = lm(formula = SoA~mag, data = SczSubset[SczSubset$domain_name == "Spatial",])
    
    # For each taken out participant (using indexes to perform same actions for HC and SCZ)
    for (LeaveOutIndex in 1:LeaveOutCount)
    {
      # Difference between coefficient of HC subject and coefficient of HC group (temporal)
      HcHcTemporalDiff = abs(sub_lm[sub_lm$sub_name == HcIds[LeaveOutIndex], "Temporal_lm_slope"]-HcTemporalModel$coefficients[2])
      
      # Difference between coefficient of HC subject and coefficient of SCZ group (temporal)
      HcSczTemporalDiff = abs(sub_lm[sub_lm$sub_name == HcIds[LeaveOutIndex], "Temporal_lm_slope"]-SczTemporalModel$coefficients[2])
      
      # Difference between coefficient of HC subject and coefficient of HC group (spatial)
      HcHcSpatialDiff = abs(sub_lm[sub_lm$sub_name == HcIds[LeaveOutIndex], "Spatial_lm_slope"]-HcSpatialModel$coefficients[2])
      
      # Difference between coefficient of HC subject and coefficient of SCZ group (spatial)
      HcSczSpatialDiff = abs(sub_lm[sub_lm$sub_name == HcIds[LeaveOutIndex], "Spatial_lm_slope"]-SczSpatialModel$coefficients[2])
      
      # Difference between coefficient of SCZ subject and coefficient of HC group (temporal)
      SczHcTemporalDiff = abs(sub_lm[sub_lm$sub_name == SczIds[LeaveOutIndex], "Temporal_lm_slope"]-HcTemporalModel$coefficients[2])
      
      # Difference between coefficient of SCZ subject and coefficient of SCZ group (temporal)
      SczSczTemporalDiff = abs(sub_lm[sub_lm$sub_name == SczIds[LeaveOutIndex], "Temporal_lm_slope"]-SczTemporalModel$coefficients[2])
      
      # Difference between coefficient of SCZ subject and coefficient of HC group (spatial)
      SczHcSpatialDiff = abs(sub_lm[sub_lm$sub_name == SczIds[LeaveOutIndex], "Spatial_lm_slope"]-HcSpatialModel$coefficients[2])
      
      # Difference between coefficient of SCZ subject and coefficient of SCZ group (spatial)
      SczSczSpatialDiff = abs(sub_lm[sub_lm$sub_name == SczIds[LeaveOutIndex], "Spatial_lm_slope"]-SczSpatialModel$coefficients[2])
      
      # We should distinguish SCZ from HC, thus:
      # Identifying SCZ as SCZ = TP (Hit)               : CorrectScz
      # Identifying SCZ as HC  = FN (Miss)              : IncorrectScz
      # Identifying HC as SCZ  = FP (False alarm)       : IncorrectHc
      # Identifying HC as HC   = TN (Correct rejection) : CorrectHc
      # To calculate accuracy: (TP+TN)/(TP+TN+FN+FP)
      
      # Sum of HC subject's coefficients differences is closer to HC group's sum of differneces than SCZ group's sum of differneces
      if ((HcHcTemporalDiff + HcHcSpatialDiff) < (HcSczTemporalDiff + HcSczSpatialDiff))
      {
        CorrectHc %+=% 1
        sub_vector[sub_vector$sub_name == HcIds[LeaveOutIndex], "Correct"] %+=% 1
        iteration_sub_vector[iteration_sub_vector$sub_name == HcIds[LeaveOutIndex], "Correct"] %+=% 1
      }
      # Sum of HC subject's coefficients differences is closer to SCZ group's sum of differneces than HC group's sum of differneces
      else
      {
        IncorrectHc %+=% 1
        sub_vector[sub_vector$sub_name == HcIds[LeaveOutIndex], "Incorrect"] %+=% 1
        iteration_sub_vector[iteration_sub_vector$sub_name == HcIds[LeaveOutIndex], "Incorrect"] %+=% 1
      }
      
      # Sum of SCZ subject's coefficients differences is closer to SCZ group's sum of differneces than HC group's sum of differneces
      if ((SczSczTemporalDiff + SczSczSpatialDiff) < (SczHcTemporalDiff + SczHcSpatialDiff))
      {
        CorrectScz %+=% 1
        sub_vector[sub_vector$sub_name == SczIds[LeaveOutIndex], "Correct"] %+=% 1
        iteration_sub_vector[iteration_sub_vector$sub_name == SczIds[LeaveOutIndex], "Correct"] %+=% 1
      }
      # Sum of SCZ subject's coefficients differences is closer to HC group's sum of differneces than SCZ group's sum of differneces
      else
      {
        IncorrectScz %+=% 1
        sub_vector[sub_vector$sub_name == SczIds[LeaveOutIndex], "Incorrect"] %+=% 1
        iteration_sub_vector[iteration_sub_vector$sub_name == SczIds[LeaveOutIndex], "Incorrect"] %+=% 1
      }
    }
    
    iterationAccuracy[IterationIndex] = sum(iteration_sub_vector$Correct) / (sum(iteration_sub_vector$Correct) + sum(iteration_sub_vector$Incorrect))
  }
  
  print(sprintf("Correct HC: %.3f, Correct PSY: %.3f, Classifier accuracy: %.3f [%.3f,%.3f], Current Time: %s", CorrectHc/(CorrectHc + IncorrectHc), CorrectScz/(CorrectScz + IncorrectScz), (CorrectScz+CorrectHc)/(CorrectScz+CorrectHc+IncorrectScz+IncorrectHc), CI(iterationAccuracy)[3], CI(iterationAccuracy)[1], format(Sys.time(), "%X")))
  
  sub_vector$CI = sprintf("[%.3f,%.3f]", CI(iterationAccuracy)[3], CI(iterationAccuracy)[1])
  
  # Returning a subject-by-subject correct/incorrect count
  return(list(sub_vector, iterationAccuracy))
}

SDTGroupCorrelations <- function()
{
  # Loading data
  exp2 <- read.csv("../Data/All_Participants_Trial_By_Trial_Data.csv", TRUE)
  
  # Confidence scores direction is inverse in the file
  exp2$conf = -1*exp2$conf
  
  exp2$group<-ifelse(str_detect(exp2$sub_name, 'Amit_HC'),"Healthy","Psychosis")
  
  exp2$sub_name = str_c(exp2$group, extract_numeric(exp2$sub_name), sep = " ")
  
  exp2$domain_name<-ifelse(exp2$domain == 1, "Temporal","Spatial")
  
  IndividualSdt = data.frame(matrix(ncol = 9, nrow = length(unique(IndividualSoA$sub_name))))
  colomnNames <- c("sub_name","TemporalHits", "TemporalMiss", "TemporalFalseAlarms", "TemporalCorrectRej", "SpatialHits","SpatialMiss", "SpatialFalseAlarms", "SpatialCorrectRej")
  colnames(IndividualSdt) <- colomnNames
  IndividualSdt$sub_name = as.factor(unique(exp2$sub_name))
  
  foreach(sub=levels(IndividualSdt$sub_name)) %do% 
    {
      TempHits = sum(exp2[exp2$sub_name == sub & exp2$domain_name == "Temporal" & exp2$mag == 0, "response"])
      TempMiss = length(exp2[exp2$sub_name == sub & exp2$domain_name == "Temporal" & exp2$mag == 0, "response"]) - TempHits
      TempFalseAlarms = sum(exp2[exp2$sub_name == sub & exp2$domain_name == "Temporal" & exp2$mag != 0, "response"])
      TempCorrectRejections = length(exp2[exp2$sub_name == sub & exp2$domain_name == "Temporal" & exp2$mag != 0, "response"]) - TempFalseAlarms
      
      SpatHits = sum(exp2[exp2$sub_name == sub & exp2$domain_name == "Spatial" & exp2$mag == 0, "response"])
      SpatMiss = length(exp2[exp2$sub_name == sub & exp2$domain_name == "Spatial" & exp2$mag == 0, "response"]) - SpatHits
      SpatFalseAlarms = sum(exp2[exp2$sub_name == sub & exp2$domain_name == "Spatial" & exp2$mag != 0, "response"])
      SpatCorrectRejections = length(exp2[exp2$sub_name == sub & exp2$domain_name == "Spatial" & exp2$mag != 0, "response"]) - SpatFalseAlarms
      
      IndividualSdt[IndividualSdt$sub_name == sub, 2:9] = c(TempHits, TempMiss, TempFalseAlarms, TempCorrectRejections, SpatHits, SpatMiss, SpatFalseAlarms, SpatCorrectRejections)
    }
  
  IndividualSdt[is.na(IndividualSdt)] <- 0
  
  SdtMeasuresTemp = dprime(n_hit = IndividualSdt$TemporalHits,
                           n_fa = IndividualSdt$TemporalFalseAlarms,
                           n_miss = IndividualSdt$TemporalMiss,
                           n_cr = IndividualSdt$TemporalCorrectRej)
  
  
  SdtMeasuresSpat = dprime(n_hit = IndividualSdt$SpatialHits, 
                           n_fa = IndividualSdt$SpatialFalseAlarms, 
                           n_miss = IndividualSdt$SpatialMiss, 
                           n_cr = IndividualSdt$SpatialCorrectRej)
  
  
  IndividualSdt$TempDprime = SdtMeasuresTemp$dprime
  IndividualSdt$TempC = SdtMeasuresTemp$c
  IndividualSdt$SpatDprime = SdtMeasuresSpat$dprime
  IndividualSdt$SpatC = SdtMeasuresSpat$c
  
  # To be able to correlate
  IndividualSdt[IndividualSdt == Inf] = 1000000
  IndividualSdt[IndividualSdt == -Inf] = -1000000
  
  HealthyCorr = corstars(as.matrix(dplyr::select(IndividualSdt[str_detect(IndividualSdt$sub_name, "Healthy"),], "TempDprime","TempC","SpatDprime","SpatC")), method = "pearson")
  
  PsychosisCorr = corstars(as.matrix(dplyr::select(IndividualSdt[str_detect(IndividualSdt$sub_name, "Psychosis"),], "TempDprime","TempC","SpatDprime","SpatC")),  method = "pearson")
  
  return(c(HealthyCorr,PsychosisCorr))
  
}

SPQB_Correlation <- function(HealthyDS, WhatToCorrelate)
{
  SPQB_DS <- read.csv("../Data/Healthy_SPQB_tbl.csv", TRUE)
  
  CorrColumnIndexes = match(WhatToCorrelate,colnames(HealthyDS))
  
  corr <- data.frame(matrix(ncol = 3, nrow = 4*length(CorrColumnIndexes)))
  cols <- c("variables", "corr", "pVal")
  colnames(corr) <- cols
  
  iter = 1
  
  for (varIndex in CorrColumnIndexes) {
    a = correlate(HealthyDS[,varIndex], SPQB_DS$SPQB_Cognitive)
    b = correlate(HealthyDS[,varIndex], SPQB_DS$SPQB_Interpersonal)
    c = correlate(HealthyDS[,varIndex], SPQB_DS$SPQB_Organization)
    d = correlate(HealthyDS[,varIndex], SPQB_DS$SPQB_Total)
    
    corr[(iter-1)*4+1,] = c(str_c("SPQ_Cog", colnames(HealthyDS)[varIndex], sep = " - "), a$correlation, a$p.value)
    corr[(iter-1)*4+2,] = c(str_c("SPQ_IntPer", colnames(HealthyDS)[varIndex], sep = " - "), b$correlation, b$p.value)
    corr[(iter-1)*4+3,] = c(str_c("SPQ_Org", colnames(HealthyDS)[varIndex], sep = " - "), c$correlation, c$p.value)
    corr[(iter-1)*4+4,] = c(str_c("SPQ_Total", colnames(HealthyDS)[varIndex], sep = " - "), d$correlation, d$p.value)
    
    iter <- iter + 1
  }
  
  return (corr)
}

numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
}

# # function for adding column of SDT category of trails -----
find_SDT<- function(mag,acc){
  if (is.na(mag)|| is.na(acc)){
    val<-"NaN"
  } else if (mag==0 && acc==1){
    val<-"hit"
  } else if (mag==0 && acc==0) {
    val<-"miss"
  } else if (mag!=0 && acc==1){
    val<-"cr"
  } else if (mag!=0 && acc==0){
    val<-"fa"
  }
}

get_plots_theme <- function(){
  val <- theme(axis.text.x = element_text(size = 13,angle = 0, vjust = 0, hjust=0.5),
               axis.text.y = element_text  (size = 13, vjust = 0.4),
               axis.title=element_text(size=14),
               panel.grid.major.y = element_line(colour = "white"),
               panel.grid.minor.y = element_line(colour = "gray"),
               panel.grid.minor.x = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.major = element_blank(),
               panel.background = element_rect(fill = "white"),
               axis.line = element_line(colour = "white"),
               plot.background=element_rect(fill="transparent",colour=NA),
               legend.key = element_rect(fill = "transparent", colour = "transparent"),
               plot.title = element_text(hjust = 0.5, size = 16),
               legend.text = element_text(size = 12),
               legend.title = element_text(size = 14))
}

corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex"), bonferroni_correcion_factor = 1)
{
  #Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P * bonferroni_correcion_factor# Matrix of p-value
  
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
  
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex") 
  }
}