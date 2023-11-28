adonis.pairwise <- function(formula, data, compare){
  
  data <- as.data.frame(data)
  
  distance <-  get(as.character(formula[[2]]))
  terms <- attr(terms.formula(formula), "term.labels")
  nTerms = length(terms)
  interTerms = grep(':', terms)
  dataIn <- data[, terms[!grepl(':', terms)], drop = F ]
  
  if(nTerms == 1){
    message('one way design')
    factorIn <- dataIn[, compare]
    if (!is.factor(factorIn)) {
      factorIn <- factor(factorIn)
    }
  }
  
  if(nTerms > 1 & length(interTerms) == 0){
    message('mutlilevel design with no interactions')
    factorIn <- dataIn[, compare]
    otherFactors <- terms[!grepl(compare, terms)]
    for(i in otherFactors)
      assign(i, dataIn[, otherFactors])
  }
  if(nTerms > 1 & length(interTerms) > 0){
    message('mutlilevel design with interactions')
    factorIn <- interaction(dataIn)
  }
  
  # Make a contrast matrix
  pairwiseRows = t(combn(levels(factorIn), 2))
  
  contrastMatrix = matrix(0, ncol = nrow(pairwiseRows), nrow = length(levels(factorIn)))
  row.names(contrastMatrix) <- levels(factorIn)
  colnames(contrastMatrix) <- apply(pairwiseRows, 1, paste, collapse = ' vs ')
  
  for(i in 1:nrow(pairwiseRows)){
    contrastMatrix[ pairwiseRows[i, 1], i] <- 1
    contrastMatrix[ pairwiseRows[i, 2], i] <- -1
  }
  
  
  # Set observations to contrast matrix values (1 or -1)
  pairwiseMatrix = 
    sapply(colnames(contrastMatrix), function(x){
      Type.temp <- factorIn
      contrasts(Type.temp) <- contrastMatrix[,x]
      Model.matrix.temp = model.matrix( ~ Type.temp)[,2]
    })
  
  pairwiseMatrix = data.frame(pairwiseMatrix)
  
  # Need to add enough factors for the residual DF to be correct,
  # So always include the same as the original anova 
  # (will repeat comparisons but they will be removed)
  Nlevels = nlevels(factorIn) - 1
  
  # Numer of comparisons to be made
  numPairwiseComps = 1:ncol(pairwiseMatrix)
  
  if(length(numPairwiseComps) == 1) stop ('No pairwise comparisons needed for only 2 groups')
  
  # Empty vectors for results
  PWadonis = NULL
  PWRes = NULL
  
  # Do the ADONIS
  for(i in numPairwiseComps){
    
    # Set factors to include
    facColsIn = c(i,numPairwiseComps[-i][1:Nlevels])
    
    # Do adonis
    if(nTerms == 1){
      # One-way
      formulaIn <- as.formula(paste(formula[[2]], "~ ."))     
      adonis.temp = adonis2(formulaIn, data = pairwiseMatrix[,facColsIn, drop = F], permutations = 999)
    }
    if(nTerms > 1 & length(interTerms) == 0){
      formulaIn <- as.formula(sprintf(paste(formula[[2]], "~ %s"), paste(paste(otherFactors, collapse="+"), '+ .')))
      adonis.temp = adonis2(formulaIn, data = pairwiseMatrix[,facColsIn, drop = F], permutations = 999)
    }
    if(nTerms > 1 & length(interTerms) > 0){
      formulaIn <- as.formula(paste(formula[[2]], "~ ."))     
      adonis.temp = adonis2(formulaIn, data = pairwiseMatrix[,facColsIn, drop = F], permutations = 999)
    }
  
  # Record Residuals
  PWRes = rbind(PWRes, adonis.temp["Residual",])
  
  # Record results
  if(nTerms == 1){
    PWadonis = rbind(PWadonis, adonis.temp[1,])
  }
  if(nTerms > 1  & length(interTerms) == 0){
    # Record block, then the results
    if(i == 1)
      PWadonis = rbind(PWadonis, adonis.temp[1,]) # block
    if(i != ncol(pairwiseMatrix))
      PWadonis = rbind(PWadonis, adonis.temp[2,]) # comaprison
  }
  if(nTerms > 1 & length(interTerms) > 0){
    PWadonis = rbind(PWadonis, adonis.temp[1,])
  }
}

as.matrix(rbind(PWadonis, adonis.temp["Residual",]))
#return(list(pairwiseTable = PWadonis, residuals = PWRes, pairwiseMatrix = pairwiseMatrix))

}
