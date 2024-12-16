
# arguments

# data: a date frame or equivalent, all cols factors
# abThresh: if corr(a,b) < this value, then delete ab arc
# abcThresh: if corr(a,b|c) < this value, then delete ab arc

pcThresh <- function(data,abThresh,abcThresh,
   myCorr='condinformation', permuteCols=FALSE)
{
   if (myCorr == 'condinformation') myCorr <- condinformation

   n <- ncol(data)

   if (permuteCols) data <- data[,sample(1:n,n,replace=FALSE)]
   adj <- matrix(1,nrow=n,ncol=n)
   adj[row(adj) > col(adj)] <- 0

   # delete unconditionally weak links
   for (i in 1:n) 
      for (j in setdiff(1:n,i))  {
         tmp <- myCorr(data[,i],data[,j])
         if (tmp < abThresh) adj[i,j] <- 0
      }

   # delete conditionally weak link
   for (k in 1:n) {
      for (i in setdiff(1:n,k)) {
         for (j in setdiff(1:n,c(k,i))) {
            if (i < j && adj[k,i] && adj[k,j] && adj[i,j]) {
               tmp <- myCorr(data[,i],data[,j],data[,k])
               if (tmp < abThresh) adj[i,j] <- 0
            }
         }
      }
   }

   diag(adj) <- 0
   colnames(adj) <- names(data)

   adj

}

