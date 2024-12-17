
# simple variant of pc algorithm, but with the big distinction that
# p-values and testing aren't used; users just set thresholds

# arguments

# data: a date frame or equivalent, all cols factors
# abThresh: if corr(a,b) < this value, then delete ab arc
# abcThresh: if corr(a,b|c) < this value, then delete ab arc

pcThresh <- function(data,abThresh,abcThresh,
   myCorrABC='condinformation', 
   myCorrAB=myCorrABC,
   outputVars=NULL,
   inputVars=NULL,
   permuteCols=FALSE)
{
   if (myCorrAB == 'condinformation') myCorrAB <- condinformation
   if (myCorrABC == 'condinformation') myCorrABC <- condinformation

   n <- ncol(data)

   varNames <- names(data)
   if (permuteCols) {
      newColNums <- sample(1:n,n,replace=FALSE)
      data <- data[,newColNums]
      newColNames <- names(data)
      # the variable varNames[i] is now in 
      # column newColNums[i]
   } else {
     newColNums <- 1:n
     newColNames <- varNames
   }



   adj <- matrix(1,nrow=n,ncol=n)
   rownames(adj) <- newColNames
   colnames(adj) <- rownames(adj)

   # e.g. want 'wageinc' to only have arrows into it, not out
   if (!is.null(outputVars)) {
      # these nodes should not to anything 
      adj[outputVars,] <- 0
   }

   # want "causation" to make sense, e.g. gender may cause occupation
   # but not vice versa
   if (!is.null(outputVars)) {
      # these nodes should not to anything 
      adj[,inputVars] <- 0
   }

   # want to have no two-way links
   for (i in 1:(n-1)) {
      if (!(i %in% outputVars)) {
         for (j in (i+1):n) {
            if (runif(1) < 0.5) adj[i,j] <- 0
            else adj[j,i] <- 0
         }
      }
   }

   # delete unconditionally weak links
   for (i in 1:n) 
      for (j in setdiff(1:n,i))  {
         tmp <- myCorrAB(data[,i],data[,j])
         if (tmp < abThresh) adj[i,j] <- 0
      }

   # delete conditionally weak link
   for (k in 1:n) {
      for (i in setdiff(1:n,k)) {
         for (j in setdiff(1:n,c(k,i))) {
            if (i < j && adj[k,i] && adj[k,j] && adj[i,j]) {
               tmp <- myCorrABC(data[,i],data[,j],data[,k])
               if (tmp < abThresh) adj[i,j] <- 0
            }
         }
      }
   }

   diag(adj) <- 0
   colnames(adj) <- names(data)

   adj

}

# forms a list l, such that l[[[i]] consists of the indices in an R
# factor f for which f == i, the i-th level of f
rowsByValue <- function(f)
{
   lapply(levels(f),function(lvl) which(f == lvl))
}

