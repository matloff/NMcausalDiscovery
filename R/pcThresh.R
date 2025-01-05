
# simple variant of pc algorithm, but with the big distinction that
# p-values and testing aren't used; users just set thresholds

# arguments

# data: a date frame or equivalent, all cols factors; if continuous
#   numeric, call regtools::discretize, e.g. discretize(x,nbins=5)
# abThresh: if corr(a,b) < this value, then delete ab arc; generally
#    small values work better, say 0.01
# abcThresh: if corr(a,b|c) < this value, then delete ab arc; again,
#    small values generally better
# outputVars: should have incoming arcs, no outgoing
# inputVars: should have outgoing arcs, no incoming, e.g.
#    avoid occupation "causing" gender
# permuteCols: one criticism of the pc alg is dependency on col order,
# autoPlot: if TRUE, plot without first saving

# needs pkgs 

#    infotheory, for default corr ftn and 'discretize'
#    igraph, for graphing
#    gtools, for writing macros
#    LogLinFit, for log-linear model (https://github.com/matloff/LogLinFit) 

pcThresh <- function(data,abThresh,abcThresh,
   outputVars=NULL,
   inputVars=NULL,
   permuteCols=FALSE,
   autoPlot=TRUE,
   logLinDegree=NULL)
{
   n <- ncol(data)
   varNames <- names(data)

   # the default and like other corr functions of interest here will
   # need all cols as factors; call discretize

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

   # use log-linear model?
   if (!is.null(logLinDegree)) {
      require(LogLinFit)
      useLogLin <- TRUE
      llOut <- llFit(data,logLinDegree)
      fitDF <- as.data.table(llOut$ary)
      # convert from cell frequences to cell probabilities
      tmp <- fitDF$value
      fitDF$value <- 1/sum(tmp) * tmp
   } else useLogLin <- FALSE

   # adjacency matrix
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
   if (!is.null(inputVars)) {
      # these nodes should not lead to anything 
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
         ### tmp <- myCorrAB(data[,i],data[,j])
         tmp <- getMutInf(i,j)
         if (tmp < abThresh) adj[i,j] <- 0
      }

   # delete conditionally weak links
   for (k in 1:n) {  # possible confounder in forks
      for (i in setdiff(1:n,k)) {
         for (j in setdiff(1:n,c(k,i))) {
            if (adj[k,i] && adj[k,j] && adj[i,j]) {
               ### tmp <- myCorrABC(data[,i],data[,j],data[,k])
               tmp <- getCondMutInf(i,j,k)
               if (tmp < abThresh) adj[i,j] <- 0
            }
         }
      }
   }

   diag(adj) <- 0

   if (autoPlot) {
      require(igraph)
      g <- graph_from_adjacency_matrix(adj,mode='directed')
      plot(g)
   }

   adj

}

# discretizes all numeric/integer columns of d, in-place
realToDiscreteFactor <- defmacro(d,
   expr={
      dscz <- infotheo::discretize
      cols <- 1:ncol(d)
      for (i in cols) {
         di <- d[,i]
         if (inherits(di,'integer')) di <- as.numeric(di)
         if (inherits(di,'numeric')) {
            d[,i] <- dscz(di,nbins=5)
         }
      }
   }
)

get1DCellProbs <- defmacro(dummy, # a macro needs an argument
   expr={
   x
   }
)

getMutInf <- defmacro(i,j,
   expr={
      tmp <- if (!useLogLin)mutinformation(data[,i],data[,j]) else
         theoretMutInfo(i,j)
   }
)

theoretMutInfo <- defmacro(i,j,
   expr={
     x
   }
)



