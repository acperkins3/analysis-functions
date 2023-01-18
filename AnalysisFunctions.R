# Analysis functions that are specific to my field trials
# For some reason, asreml won't take formulas, so we need this eval(parse form

'%!in%' <- function(x,y)!('%in%'(x,y))

# na.method needs to be included for the models with the autoregressive variance structures


GetStandCountpValue <- function(data, Trait) {
  if (!is.na(data[1, "StandCount"])){
    SubblockTerm <- ifelse(data[1, "Germplasm"] == "GEMN SSD X LH244", " + Replicate:SubBlock", "")
    model <- eval(parse(text = paste0("asreml(fixed = ", Trait, " ~ Pedigree + StandCount, random = ~Replicate", SubblockTerm, ", residual = ~idv(units), na.action = na.method(x = 'include', y = 'include'), data = data, trace = FALSE)")))
    pval <- wald(model)
    return(wald(model)[3,4]) # Note: This indexing is dependent on the number of fixed effects
  }
  else {
    return(NA)
  }
}

# We don't know which models didn't converge, so this function finds and removes those
# It checks the models with Pedigree as a fixed effect since the goal is to 
# use models to calculate BLUEs. The random models are just for heritability

PoliceConvergence <- function(model) {
  if (model$fixed$converge == FALSE) {
    return(NULL)
  }
  else {
    return(model)
  }
}


FitIIDModel <- function(data, Trait, StandCountP) {
  if (!is.na(StandCountP) & length(unique(data$StandCount)) > 1) { #If there are stand counts
    if (StandCountP > 0.05) { #If the effect of stand count is not significant
      data$StandCount <- NA #Make the stand counts NA so they are not used in the model (for this use of the function)
    }
  }
  StandCountTerm <- ifelse(!is.na(data[1, "StandCount"]), " + StandCount", "")
  SubblockTerm <- ifelse(data[1, "Germplasm"] == "GEMN SSD X LH244", " + Replicate:SubBlock", "")
  model <- PoliceConvergence(list(fixed = eval(parse(text = paste("asreml(fixed =", Trait, "~ Pedigree", StandCountTerm, ", random = ~Replicate", SubblockTerm, ", residual = ~idv(units), na.action = na.method(x = 'include', y = 'include'), data = data, trace = FALSE)"))),
                                  random = eval(parse(text = paste("asreml(fixed =", Trait, "~ 1",  StandCountTerm,", random = ~Pedigree + Replicate", SubblockTerm, ", residual = ~idv(units), na.action = na.method(x = 'include', y = 'include'), data = data, trace = FALSE)")))))
}

FitAR1RowModel <- function(data, Trait, StandCountP) {
  if (!is.na(StandCountP) & length(unique(data$StandCount)) > 1) { #If there are stand counts
    if (StandCountP > 0.05) { #If the effect of stand count is not significant
      data$StandCount <- NA #Make the stand counts NA so they are not used in the model (for this use of the function)
    }
  }
  StandCountTerm <- ifelse(!is.na(data[1, "StandCount"]), " + StandCount", "")
  SubblockTerm <- ifelse(data[1, "Germplasm"] == "GEMN SSD X LH244", " + Replicate:SubBlock", "")
  model <- PoliceConvergence(list(fixed = eval(parse(text = paste("asreml(fixed =", Trait, "~ Pedigree", StandCountTerm, ", random = ~Replicate", SubblockTerm, ", residual = ~ar1(Row):idv(Col), na.action = na.method(x = 'include', y = 'include'), data = data, trace = FALSE)"))),
                                  random = eval(parse(text = paste("asreml(fixed =", Trait, "~ 1",  StandCountTerm,", random = ~Pedigree + Replicate", SubblockTerm, ", residual = ~ar1(Row):idv(Col), na.action = na.method(x = 'include', y = 'include'), data = data, trace = FALSE)")))))
}

FitAR1ColModel <- function(data, Trait, StandCountP) {
  if (!is.na(StandCountP) & length(unique(data$StandCount)) > 1) { #If there are stand counts
    if (StandCountP > 0.05) { #If the effect of stand count is not significant
      data$StandCount <- NA #Make the stand counts NA so they are not used in the model (for this use of the function)
    }
  }
  StandCountTerm <- ifelse(!is.na(data[1, "StandCount"]), " + StandCount", "")
  SubblockTerm <- ifelse(data[1, "Germplasm"] == "GEMN SSD X LH244", " + Replicate:SubBlock", "")
  model <- PoliceConvergence(list(fixed = eval(parse(text = paste("asreml(fixed =", Trait, "~ Pedigree", StandCountTerm, ", random = ~Replicate", SubblockTerm, ", residual = ~idv(Row):ar1(Col), na.action = na.method(x = 'include', y = 'include'), data = data, trace = FALSE)"))),
                                  random = eval(parse(text = paste("asreml(fixed =", Trait, "~ 1",  StandCountTerm,", random = ~Pedigree + Replicate", SubblockTerm, ", residual = ~idv(Row):ar1(Col), na.action = na.method(x = 'include', y = 'include'), data = data, trace = FALSE)")))))
}

FitAR1BothModel <- function(data, Trait, StandCountP) {
  if (!is.na(StandCountP) & length(unique(data$StandCount)) > 1) { #If there are stand counts
    if (StandCountP > 0.05) { #If the effect of stand count is not significant
      data$StandCount <- NA #Make the stand counts NA so they are not used in the model (for this use of the function)
    }
  }
  StandCountTerm <- ifelse(!is.na(data[1, "StandCount"]), " + StandCount", "")
  SubblockTerm <- ifelse(data[1, "Germplasm"] == "GEMN SSD X LH244", " + Replicate:SubBlock", "")
  model <- PoliceConvergence(list(fixed = eval(parse(text = paste("asreml(fixed =", Trait, "~ Pedigree", StandCountTerm, ", random = ~Replicate", SubblockTerm, ", residual = ~ar1v(Row):ar1(Col), na.action = na.method(x = 'include', y = 'include'), data = data, trace = FALSE)"))),
                                  random = eval(parse(text = paste("asreml(fixed =", Trait, "~ 1",  StandCountTerm,", random = ~Pedigree + Replicate", SubblockTerm, ", residual = ~ar1v(Row):ar1(Col), na.action = na.method(x = 'include', y = 'include'), data = data, trace = FALSE)")))))
}


# These get the correlations in the row and column directions from models that include a spatial correction

GetRowCor <- function(BestModel) {
  round(summary(BestModel$fixed)$varcomp["Row:Col!Row!cor", "component"], 2)
}

GetColCor <- function(BestModel) {
  round(summary(BestModel$fixed)$varcomp["Row:Col!Col!cor", "component"], 2)
}



# The asreml row and column variance functions don't work when not all the row x column values
# are included. Hence, expand.grid adds the missing ones
# The data need to be sorted by row and then column, but I think that is done
# by the merge function

ExpandData <- function(data) {
  dataexpanded <- expand.grid(Row = 1:max(as.numeric(as.character(data$Row)), na.rm=TRUE), Col = 1:max(as.numeric(as.character(data$Col)), na.rm=TRUE)) %>%
    mutate(Row = as.factor(Row)) %>%
    mutate(Col = as.factor(Col))
  
  dataexpanded <- merge(dataexpanded, data, by = c("Row", "Col"), all = TRUE) %>%
    mutate(Row = as.factor(Row)) %>%
    mutate(Col = as.factor(Col))
}

# BIC is from the model with genotype as a fixed effect since the goal is to get BLUEs for second stage

GetBIC <- function(model) {
  if (is.null(model)) {
    return(NA)
  }
  else {
    return(summary(model$fixed)$bic[1])
  }
}

# This is to find the best model (based on BIC) and return the best model object (in a new column)
#Lower BIC is better, so we fill in the missing stuff with a really large value
#Might there be a tie?

GetBestModel <- function(IIDModel, IIDBIC, AR1RowModel, AR1RowBIC, AR1ColModel, AR1ColBIC, AR1BothModel, AR1BothBIC) {
  if(is.na(IIDBIC)) {IIDBIC <- 1000000}
  if(is.na(AR1RowBIC)) {AR1RowBIC <- 1000000}
  if(is.na(AR1ColBIC)) {AR1ColBIC <- 1000000}
  if(is.na(AR1BothBIC)) {AR1BothBIC <- 1000000}
  if(IIDBIC < AR1RowBIC & IIDBIC < AR1ColBIC & IIDBIC < AR1BothBIC) {return(IIDModel)}
  if(AR1RowBIC < IIDBIC & AR1RowBIC < AR1ColBIC & AR1RowBIC < AR1BothBIC) {return(AR1RowModel)}
  if(AR1ColBIC < IIDBIC & AR1ColBIC < AR1RowBIC & AR1ColBIC < AR1BothBIC) {return(AR1ColModel)}
  if(AR1BothBIC < IIDBIC & AR1BothBIC < AR1RowBIC & AR1BothBIC < AR1ColBIC) {return(AR1BothModel)}
  
}

GetBLUEs <- function(model, data) {
  if (!is.null(model)) {
    model <- eval(model$fixed$call)
    predictions <- predict(model, classify = "Pedigree", vcov = TRUE)
    
    pvals <- as.data.frame(predictions$pvals)
    
    vcov <- as.matrix(predictions$vcov) #This is taken from the ASRtriala source code
    sel <- matrix(1, ncol = 1, nrow = length(pvals$predicted.value))
    sel[is.na(pvals$predicted.value), ] <- 0
    vcov <- vcov[sel == 1, sel == 1]
    pvals$Weight[sel == 1] <- diag(solve(vcov))
    
    for (row in 1:nrow(pvals)) { #This is clunky and should be replaced someday
      temp <- data %>% filter(!is.na(ThisTrait)) %>%
        mutate(Pedigree = as.character(Pedigree)) %>%
        filter(Pedigree == pvals[row, "Pedigree"])
      pvals[row, "SampleSize"] <- nrow(temp)
    }
    
    return(pvals)
  }
  else {
    return(NULL)
  }
}

# The predict function returns an error sometimes because ASReml-R runs out of memory
# so the use of `safely` here just makes a note of that rather than exiting the map
# function

GetBLUEsSafely <- safely(GetBLUEs)


GetBLUPs <- function(model, data) {
  if (!is.null(model)) {
    model <- eval(model$random$call)
    return(as.data.frame(predict(model, classify = "Pedigree")$pvals))
  }
  else {
    return(NULL)
  }
}

GetBLUPsSafely <- safely(GetBLUPs)

AddTraitColumnToData <- function(data, trait) {
  data[,"ThisTrait"] <- data[,trait]
  return(data)
}


GetGeneticVariance <- function(model) {
  summary(model$random)$varcomp["Pedigree","component"]
}

# SED Function

GetSED <- function(model, data) {
  if (!is.null(model$random)) {
    model <- eval(model$random$call)
    return(as.numeric(predict(model, classify = "Pedigree")$avsed))
  }
  else {
    return(NULL)
  }
}

GetSEDSafely <- safely(GetSED)

# New H2 Function

GetCullisH2SED <- function(SED, GeneticVariance) {
  if (is.null(SED$result)) {
    return(NULL)
  }
  else {
    SED <- as.numeric(SED$result)
    Cullis <- 1 - (SED**2 / (2*GeneticVariance))
    return(round(Cullis, 3))
  }
}



