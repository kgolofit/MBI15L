test <- function() {
  require(Biostrings)
  
  seq1 = sample(DNA_ALPHABET[1:4], size=30, replace=TRUE)
  seq2 = sample(DNA_ALPHABET[1:4], size=10, replace=TRUE)
  
  doGotoh(seq1, seq2, local=FALSE)
  print('=======================================')
  pairwiseAlignment(paste(seq1, collapse=" "), paste(seq2, collapse=" "), gapOpening=-12, gapExtension=-10)
  
#   doGotoh(c('A', 'A', 'A','G', 'G', 'T', 'T'), c('A', 'A', 'A', 'T', 'T'))
  doAltschulErickson(seq1, seq2, 10, 12, 10)
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
## seq1 - sequence 1
## seq2 - sequence 2
## v, u - affine gap cost coeff.
## dist - distance between different chars
doAltschulErickson <- function(seq1, seq2, v, u, dist) {
  M = length(seq1)
  N = length(seq2)
  
  print(M)
  print(N)
  
  ## Step {1}
  P = matrix(nrow=M, ncol=N)
  Q = matrix(nrow=M, ncol=N)
  R = matrix(nrow=M, ncol=N)
  a = matrix(0, nrow=M+1, ncol=N+1)
  b = matrix(0, nrow=M+1, ncol=N+1)
  c = matrix(0, nrow=M+1, ncol=N+1)
  d = matrix(0, nrow=M+1, ncol=N+1)
  e = matrix(0, nrow=M+1, ncol=N+1)
  f = matrix(0, nrow=M+1, ncol=N+1)
  g = matrix(0, nrow=M+1, ncol=N+1)
  
  for(j in 1:N) {
    P[1,j] = Inf
    R[1,j] = v + j*u
  }
  for(i in 1:M) {
    Q[i,1] = Inf
    R[i,1] = v + i*u
  }
  R[1,1] = 0
  
  c[M+1,N+1] = 1
  
  for(i in 2:M) {
    for(j in 2:N) {
      if (seq1[i] == seq2[j]) {
        cost = 0
      }
      else {
        cost = dist
      }
      ## Step {2}
      P[i,j] = u + min(P[i-1,j], R[i-1,j]+v)
      
      ## Step {3}
      if(P[i,j] == P[i-1,j]+u) {
        d[i-1,j] = 1
      }
      if(P[i,j] == R[i-1,j]+v+u) {
        e[i-1,j] = 1
      }
      
      ## Step {4}
      Q[i,j] = u + min(Q[i,j-1], R[i,j-1]+v)
      
      ## Step {5}
      if(Q[i,j] == Q[i,j-1]+u) {
        f[i,j-1] = 1
      }
      if(Q[i,j] == R[i,j-1]+v+u) {
        g[i,j-1] = 1
      }
      
      ## Step {6}
      R[i,j] = min(P[i,j], Q[i,j], R[i-1,j-1]+cost)
      
      ## Step {7}
      if(R[i,j] == P[i,j]) {
        a[i,j] = 1
      }
      if(R[i,j] == Q[i,j]) {
        b[i,j] = 1
      }
      if(R[i,j] == R[i-1,j-1]+cost) {
        c[i,j] = 1
      }
    }
  }
  
  for(i in seq(M, 2, by=-1)) {
    for(j in seq(N, 2, by=-1)) {
      ## Step {8}
      if((a[i+1,j]==0 || e[i,j]==0) && (b[i,j+1]==0 || g[i,j]==0) && c[i+1,j+1]==0) {
        a[i,j] = 0
        b[i,j] = 0
        c[i,j] = 0
      }
      
      ## Step {9}
      if(a[i+1,j]==0 && b[i,j+1]==0 && c[i+1,j+1]==0) {
        next
      }
      
      ## Step {10}
      if(a[i+1,j]==1 && d[i,j]==1) {
        d[i+1,j] = 1-e[i,j]
        e[i,j] = 1-a[i,j]
        a[i,j] = 1
      }
      else {
        d[i+1,j] = 0
        e[i,j] = 0
      }
      
      ## Step {11}
      if(b[i,j+1]==1 && f[i,j]==1) {
        f[i,j+1] = 1-g[i,j]
        g[i,j] = 1-b[i,j]
        b[i,j] = 1
      }
      else {
        f[i,j+1] = 0
        g[i,j] = 0
      }
    }
  }
  
  print(R)
  #print(a)
  #print(b)
  #print(c)
  #print(d)
  #print(e)
  #print(f)
  #print(g)
}