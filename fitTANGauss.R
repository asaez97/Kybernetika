# Conditional Gaussian TAN
# Chow-Liu Algorithm--------------------
## Compute TAN structure----------------------------
fit_tan_structure = function(mutualInfoCond,target,root=1,pred_disc){
  # Compute maximun spanning tree
  
  # Number of nodes
  n <- nrow(mutualInfoCond)
  if(is.null(n)){
    n=1
  }
  # Inicializate
  selected_nodes <- rep(FALSE, n)
  names(selected_nodes) = row.names(mutualInfoCond)
  # root of maximal tree
  selected_nodes[root] <- TRUE  
  # Adjacency matrix of tree
  maximal_tree <- matrix(0, n, n,dimnames = dimnames(mutualInfoCond))  
  
  # Compute maximal spanning tree
  i = 1
  for (i in 1:(n-1)) {
    max_weight <- -Inf
    u <- -1#row
    v <- -1# Column
    if(all(selected_nodes[pred_disc])){
      # All discrete nodes are selected
      option_nodes = names(which(!selected_nodes))
    }else{
      # There is not discrete nodes selected
      option_nodes = names(which(!selected_nodes[pred_disc]))
    }
    # Find the node which is not selected with better value
    for (j in names(which(selected_nodes))) {
      for (k in option_nodes) {
        if (mutualInfoCond[j, k] > max_weight) {
          max_weight <- mutualInfoCond[j, k]
          u <- j
          v <- k
        }
      }
    }
    
    # Add the node to the tree
    if (u != -1 && v != -1) {
      maximal_tree[u, v] <- 1
      selected_nodes[v] <- TRUE
    }
  }
  # variable names
  varTAN <- c(row.names(maximal_tree),target)
  # Adding the class to the tree as father
  adj_mat <- rbind(maximal_tree,1)
  adj_mat <- cbind(adj_mat,0)
  rownames(adj_mat)<- varTAN
  colnames(adj_mat)<- varTAN
  
  return(adj_mat)
}

# choose root---------------------
fit_root = function(MI,pred_disc = character()){
  
  if(length(pred_disc)!=0){
    # If there are discrete feature variabes, one of them must be the root
    MI = MI[pred_disc,pred_disc,drop = F]
  }
  root = which.max(MI)%%nrow(MI)
  root = ifelse(root==0,nrow(MI),root)# Al hacer el modulo, 0=root=nrow
  root = colnames(MI)[root]
  
  return(root)
}

# Learning CG TAN--------------------------------
fit_tan_g <- function(target, data, root = NULL,
                      mutualInfoCond = NULL, all = F){
  
  if(is.null(target)){
    stop(paste0("Variable ", taget, " must be indicated"))
  }
  
  if(!(target%in%colnames(data))){
    stop(paste0("Variable ", target, " must be a name of a variable of data"))
  }
  # Structural learning--------------------------------------
  # Compute the MI matrix
  if(is.null(mutualInfoCond)){
    mutualInfoCond = MI_tan_gauss(data,target)
  }
  # Is there a discrete variable?
  varPred = setdiff(colnames(data),target)
  pred_cont = sapply(varPred,function(variable){
    is.numeric(data[1,variable])
  })
  pred_cont = names(pred_cont)[which(pred_cont)]
  pred_disc = setdiff(varPred,pred_cont)
  # browser()
  # Compute the root
  if(length(pred_disc)>0){
    if(is.null(root)||!(root%in%pred_disc)){
     root = fit_root(mutualInfoCond,pred_disc)
    }
  }else if(is.null(root)){
    root = fit_root(mutualInfoCond,pred_disc)
  }
  
  # Compute the adjacency matrix
  adj_mat = fit_tan_structure(mutualInfoCond,target,root,pred_disc)
  dag = empty.graph(colnames(data))
  amat(dag)<-adj_mat
  bn = bn.fit(dag,data = data)
  
  if(all){# Devolver informacion mutua y distribicones tambien
    return(list(bn=bn,mutualInfoCond=mutualInfoCond))
  }else{
    return(bn)
  }
}

## Mutual information functions--------------------------------------------

cond_mi_cont = function(X,Y,target,data){
  levels_target = levels(data[[target]])
  prob_target = as.vector(prop.table(table(data[[target]])))
  
  corXY_t = sapply(levels_target, function(i){
    cor(data[data[[target]]==i,X],data[data[[target]]==i,Y])
  })
  return(-1/2*sum(prob_target*log((1-corXY_t^2))))
}

cond_mi_cont_disc = function(X,Y,target,data){
  # target or class variable distribution
  levels_target = levels(data[[target]])
  prob_target = as.vector(prop.table(table(data[[target]])))
  # Feature variable distributions
  ## Continuous
  varX_t = sapply(levels_target, function(i){
    var(data[data[[target]]==i,X])
  })
  #Discrete
  levels_Y = levels(data[[Y]])
  prob_joint = prop.table(table(data[,c(target,Y)]))
  # Compute Mutual information
  var_Y_T = c()
  i = levels_target[1]
  j = levels_Y[1]
  for(i in levels_target){
    var_Y_T = rbind(var_Y_T,
                    sapply(levels_Y,function(j){
                      var(data[(data[[target]]==i & data[[Y]]==j),X])
                    }
                    )
    )
  }
  # Variances very smalls
  var_Y_T = ifelse(var_Y_T<10^-6,10^-6,var_Y_T)
  return(1/2*(sum(prob_target*log(varX_t))-sum(prob_joint*log(var_Y_T),na.rm = T)))
}

# Mutual information of dataset------------------
MI_tan_gauss = function(data,target){
  varPred = setdiff(colnames(data),target)
  MI = matrix(0,ncol = length(data)-1,nrow = length(data)-1,
              dimnames = list(varPred,varPred))
  
  i = 1
  j = 2
  # Compute MI between each pair of features variables
  for(i in 1:(length(varPred)-1)){
    Y = varPred[i]
    Y_numeric = is.numeric(data[[Y]])
    for(j in (i+1):length(varPred)){
      X = varPred[j]
      if(is.numeric(data[[X]])){
        if(Y_numeric){
          MI[X,Y] = cond_mi_cont(X,Y,target,data)
          MI[Y,X] = MI[X,Y]
        }
      }else{
        if(Y_numeric){
          MI[X,Y] = cond_mi_cont_disc(Y,X,target,data)
          MI[Y,X] =  MI[X,Y] # CG asumptions
        }else{
          MI[X,Y] = infotheo::condinformation(data[[X]],data[[Y]],data[[target]])
          MI[Y,X] = MI[X,Y]
        }
      }
    }
  }
  return(MI)
}
