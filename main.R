test <- function() {
  require(Biostrings)
  
  seq1 = sample(DNA_ALPHABET[1:4], size=30, replace=TRUE)
  seq2 = sample(DNA_ALPHABET[1:4], size=10, replace=TRUE)
  
  doGotoh(seq1, seq2, local=FALSE)
  print('=======================================')
  pairwiseAlignment(paste(seq1, collapse=" "), paste(seq2, collapse=" "), gapOpening=-12, gapExtension=-10)
  
#   doGotoh(c('A', 'A', 'A','G', 'G', 'T', 'T'), c('A', 'A', 'A', 'T', 'T'))
    
}

#############################################
## Gotoh algorithm implementation
##
## Based on "An improved algorithm for matching biological sequences" article.
##
## Params:
## seq1 - sequence 1
## seq2 - sequence 2
## distEq - distance when chars from seq1 and seq2 are equal
## distNEq - distance when chars from seq1 and seq2 are NOT equal
## u,v - coefficiency of linear gap penality function
## local - local alignment
##
doGotoh <- function(seq1, seq2, distEq=0, distNEq=1, u=10, v=12, local=FALSE) {
  ## initialization of matrices
  dMatrix <- matrix(0, nrow=length(seq1) + 1, ncol=length(seq2) + 1)
  qMatrix <- matrix(0, nrow=length(seq1) + 1, ncol=length(seq2) + 1)
  pMatrix <- matrix(0, nrow=length(seq1) + 1, ncol=length(seq2) + 1)
  backMatrix <- matrix(0, nrow=length(seq1) + 1, ncol=length(seq2) + 1)
  
  if(!local) {
    for(i in 2:dim(dMatrix)[1]) {
      dMatrix[i,1] <- u * (i-2) + v
      pMatrix[i,1] <- u * (i-2) + v
      qMatrix[i,1] <- u * (i-2) + v
    }
    
    for(i in 2:dim(dMatrix)[2]) {
      dMatrix[1,i] <- u * (i-2) + v
      pMatrix[1,i] <- u * (i-2) + v
      qMatrix[1,i] <- u * (i-2) + v
    }
  }
  
  ## taka funkcja
  for(i in 2:dim(dMatrix)[1]) {
    for(j in 2:dim(dMatrix)[2]) {
      
      if(seq1[i-1] != seq2[j-1]) {
        dist <- distNEq
      }
      else {
        dist <- distEq
      }
      
      pMatrix[i,j] <- min(dMatrix[i-1,j] + u + v,
                          pMatrix[i-1,j] + u)
      qMatrix[i,j] <- min(dMatrix[i,j-1] + u + v,
                          qMatrix[i,j-1] + u)
      dMatrix[i,j] <- min((dMatrix[i-1, j-1] + dist),
                          pMatrix[i,j], qMatrix[i, j])
      
      
      if(dMatrix[i, j] == pMatrix[i,j]) {
        backMatrix[i,j] <- 2
      }
      else if (dMatrix[i, j] == qMatrix[i, j]) {
        backMatrix[i,j] <- 3
      }
      else if(dMatrix[i, j] == (dMatrix[i-1, j-1] + dist)) {
        backMatrix[i,j] <- 1
      }
    }
  }
  
  ## backtracking
  aString <- ""
  bString <- ""
  if(local) {
    
    i <- which.min(dMatrix[,dim(dMatrix)[2]])
    j <- which.min(dMatrix[dim(dMatrix)[1],])
    print(i)
    print(j)
  }
  else {
    i <- dim(backMatrix)[1]
    j <- dim(backMatrix)[2]
  }
  while(backMatrix[i,j] != 0) {
    if(backMatrix[i,j] == 1) {
      aString <- paste(seq1[i-1], aString)
      bString <- paste(seq2[j-1], bString)
      
      i <- i-1
      j <- j-1
    }
    else if(backMatrix[i,j] == 2) {
      aString <- paste(seq1[i-1], aString)
      bString <- paste('-', bString)
      
      i <- i-1
    }
    else if(backMatrix[i,j] == 3) {
      aString <- paste('-', aString)
      bString <- paste(seq2[j-1], bString)
      
      j <- j-1      
    }
  }
  
  print(dMatrix)
  print(backMatrix)
  
  print(seq1)
  print(seq2)
  
  print(aString)
  print(bString)
}

##############################################
## Altschul-Erickson algorithm implementation
##
doAltschulErickson <- function(seq1, seq2) {
  
}