# Funciones para crear el TAN

# Función para el algoritmo de chow-liu--------------------
## Función para el arbol maximal--------------------
prim_maximal <- function(adj_mat,root = 1) {
  # Número de nodos en el grafo
  # browser()
  n <- nrow(adj_mat)
  if(is.null(n)){
    n=1# caso 1 variable
  }
  
  # Inicialización
  selected_nodes <- rep(FALSE, n)
  names(selected_nodes) = row.names(adj_mat)
  selected_nodes[root] <- TRUE  # Seleccionamos el nodo inicial arbitrario
  maximal_tree <- matrix(0, n, n,dimnames = dimnames(adj_mat))  # Para almacenar el árbol máximo
  
  # Bucle principal
  i = 1
  for (i in 1:(n-1)) {
    max_weight <- -Inf
    u <- -1#Fila
    v <- -1# Columna
    
    # Encontrar la arista de mayor peso que conecta un nodo seleccionado con uno no seleccionado
    for (j in which(selected_nodes)) {
      for (k in which(!selected_nodes)) {
        if (adj_mat[j, k] > max_weight) {
          max_weight <- adj_mat[j, k]
          u <- j
          v <- k
        }
      }
    }
    
    # Añadir la arista al árbol maximal
    if (u != -1 && v != -1) {
      maximal_tree[u, v] <- 1
      # maximal_tree[v, u] <- max_weight
      selected_nodes[v] <- TRUE
    }
  }
  return(maximal_tree)
}

## Funcion para generar la estructura TAN----------------------------
fit_tan_structure = function(mutualInfoCond,target,root=1){
  # browser()
  # Matriz de adjacencia
  adj_mat <- prim_maximal(mutualInfoCond,root = root)
  # nombre de las variables
  varTAN <- c(row.names(adj_mat),target)
  # anadimos la clase
  adj_mat <- rbind(adj_mat,1)# Aristas del target a todos los predictores
  adj_mat <- cbind(adj_mat,0)
  rownames(adj_mat)<- varTAN
  colnames(adj_mat)<- varTAN
  
  # Calculamos el tan
  # library(bnlearn)
  # dag <- empty.graph(varTAN)
  # amat(dag)<- adj_mat
  # graphviz.plot(dag)
  return(adj_mat)
}
# Funcion para calcular la matriz de IM dada la clase---------------------------

mutual_information_tan = function(data,target,fit.args,root = 1,parallel = F){
  # Tan para regresion
  # browser()
  source("mutualInformation.R")
  # Seleccionamos las variables
  varNames = colnames(data)
  # variables predictoras
  varPred = varNames[varNames!=target]
  
  
  # Inicializamos la matriz para MI
  # Inicializamos la matriz para MI
  MI = matrix(0,nrow = length(varPred),ncol = length(varPred),
              dimnames = list(varPred,varPred))
  ff = function(i,X,target,data,distC,distX_C,distributions,fit.args){
    Y = varPred[i]
    # browser()
    if(X==Y){
      return(list(MI=0,name = Y))
    }else{
      loginfo(paste("Learning",
                    paste0(X,"|",paste0(sort(c(Y,target)),collapse = ":")),
                    collapse = " "))
      distY_C = distributions[[paste0(Y,"|",target)]]
      # print(distY_C)
      # exists("cond_mut_information_cont")
      if(distX_C$varType=="Continuous"){
        if(distY_C$varType=="Continuous"){
          result = tryCatch(cond_mut_information_cont(data = data,className = target,
                                             varNames = c(X,Y,target),
                                             fit.args = fit.args,
                                             distC = distC,distY_C = distY_C,
                                             distX_C = distX_C),
                   error = function(e) list(MI=0))
        }else{# Quien condiciona es discreto
          result = tryCatch(cond_mut_information_cont_dis(data = data,className = target,
                                                          varNames = c(X,Y,target),
                                                          fit.args = fit.args,
                                                          distC = distC,distY_C = distY_C,
                                                          distX_C = distX_C),
                            error = function(e) list(MI=0))
        }
      }else{
        # La distribucion de X es discreta
        if(distY_C$varType=="Continuous"){
          result = tryCatch(cond_mut_information_dis_cont(data = data,className = target,
                                                          varNames = c(X,Y,target),
                                                          fit.args = fit.args,
                                                          distC = distC,distY_C = distX_C,
                                                          distX_C = distY_C),
                            error = function(e) list(MI=0))
        }else{# Quien condiciona es discreto
          result = cond_mut_information_dis(data = data,className = target,
                                            varNames = c(X,Y,target),
                                            fit.args = fit.args,
                                            distC = distC,distY_C = distY_C,
                                            distX_C = distX_C)
        }
      }
      # return(list(MI=0,distributions = distY_C))
    }
    return(list(MI = result$MI,name = Y))
  }
  if(parallel){
    time = Sys.time()
    library(foreach)
    library(doParallel)
    # setup parallel backend to use many processors
    cores = detectCores()
    # cl <- makeCluster(cores[1]-1) #not to overload your computer
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)
    result = foreach(var_i = varPred,
                     .packages = "MoTBFs")%dopar%{
                       source("mutualInformation.R")
                       list(distX_cond(data = data,var_i,target,fit.args = fit.args),
                            paste0(var_i,"|",target))
                     }
    distributions = lapply(result,"[[",1)
    names(distributions)=lapply(result,"[[",2)
    distC = distVar(data[,target,drop = F],target,fit.args)
    # stopCluster(cl)
    # X = varPred[1]
    # distributions2 = list()
    # cl <- makeCluster(cores[1]-1)
    # registerDoParallel(cl)
    foreach(X = varPred,
            .packages = c("logging","MoTBFs"))%do%{
              # variable que es condicionada
              # Aprendizaje de la distribucion var_i|var_j,C
              distX_C = distributions[[paste0(X,"|",target)]]
              
              
              result = foreach(i = 1:length(varPred),
                               .packages = c("logging","MoTBFs"))%dopar%{
                                 source("mutualInformation.R")
                                 ff(i,X,target,data,distC,distX_C,
                                    distributions,fit.args)
                               }
              
              
              MIs = sapply(result, "[[", "MI")
              namesDist = sapply(result, "[[", "name")
              # distr2 = sapply(result,"[[","distributions")
              # names(distr2) = sapply(namesDist, "[", 2)
              # distributions2 = append(distr2,distributions2)
              # Padres por filas e hijos por columnas
              MI[sapply(namesDist, "[", 1),X]=MIs
              
            }
    # result=list(MI,distributions)
    # save(result,file=paste("resultados_iter_",iter,".Rdata",sep =""))
    stopCluster(cl)
    Sys.time()-time
  }else{
    # Distribucion de la variable clase y de X_i|C--------------------
    distributions = fitNB_dis(data,target,fit.args)
    # if(is.null(root)){
    #   MI = mutual_info_df_Y(data,target,fit.args,distributions)
    #   root = which.max(MI)
    # }
    X = varPred[3]
    distC = distributions[[target]]
    # nDist = length(distributions)
    # browser()
    for(X in varPred){
      # variable que es condicionada
      # Aprendizaje de la distribucion var_i|var_j,C
      distX_C = distributions[[paste0(X,"|",target)]]
      Y = varPred[2]
      for(Y in varPred){
        # iter = iter+1
        if(X==Y){
          next
        }
        # if(iter<8000){
        #   next
        # }
        
        loginfo(paste("Learning",
                      paste0(X,"|",paste0(sort(c(Y,target)),collapse = ":")),
                      collapse = " "))
        distY_C = distributions[[paste0(Y,"|",target)]]
        if(distX_C$varType=="Continuous"){
          if(distY_C$varType=="Continuous"){
            result=cond_mut_information_cont(data = data,className = target,
                                             varNames = c(X,Y,target),
                                             fit.args = fit.args,
                                             distC = distC,distY_C = distY_C,
                                             distX_C = distX_C)
          }else{# Quien condiciona es discreto
            result=cond_mut_information_cont_dis(data = data,className = target,
                                                 varNames = c(X,Y,target),
                                                 fit.args = fit.args,
                                                 distC = distC,distY_C = distY_C,
                                                 distX_C = distX_C)
          }
        }else{
          # La distribucion de X es discreta
          if(distY_C$varType=="Continuous"){
            result=cond_mut_information_dis_cont(data = data,className = target,
                                                 varNames = c(X,Y,target),
                                                 fit.args = fit.args,
                                                 distC = distC,distY_C = distX_C,
                                                 distX_C = distY_C)
          }else{# Quien condiciona es discreto
            result=cond_mut_information_dis(data = data,className = target,
                                            varNames = c(X,Y,target),
                                            fit.args = fit.args,
                                            distC = distC,distY_C = distY_C,
                                            distX_C = distX_C)
          }
        }
        # nDist=nDist+1
        # distributions[nDist]=result$distributions[1]
        # names(distributions)[nDist]=names(result$distributions[1])
        # Padres por filas e hijos por columnas
        MI[Y,X]=result$MI
        # result=list(MI,distributions)
        # save(result,file=paste("resultados_iter_",iter,".Rdata",sep =""))
      }
      save(MI,file=paste("../Experimentacion/MI.RData",sep =""))
    }
  }
  # Hacemos la media para considerar un grafo no dirigido
  MI = (MI +t(MI))/2
  if(is.null(root)){
    root = fitroot(MI)
  }
  return(list(MI=MI,root = root))
}
# Mutual Information tan discreto y CG-----------------------------
mutual_information_tan_d = function(data_new,varsNew=NULL,target,
                                    mutualInfoCond_0=NULL){
  # browser()
  # Seleccionamos las variables
  varPred = setdiff(colnames(data_new),target)
  # variables predictoras antiguas
  varsOld = setdiff(varPred,varsNew)
  
  # Distribuciones de la variable nuevas X_i|C--------------------
  
  # if(is.null(root)){
  #   MI = mutual_info_df_Y(data,target,fit.args,distributions)
  #   root = which.max(MI)
  # }
  # Inicializamos la matriz para MI
  MI_new = matrix(0,nrow = length(varPred),ncol = length(varPred),
                  dimnames = list(varPred,varPred))
  if(!is.null(mutualInfoCond_0)){
    MI_new[varsOld,varsOld]=mutualInfoCond_0
  }
  # browser()
  for(X in varsNew){
    for(Y in varPred){
      if(X==Y){
        next
      }else{
        MI = condinformation(data_new[,X],data_new[,Y],data_new[,target])
        MI_new[Y,X] = MI
        MI_new[X,Y] = MI
      }
    }
  }
  return(MI_new)
}
# funcion para calcular la raiz de un TAN---------------------
fitroot = function(MI){
  root = which.max(MI)%%nrow(MI)# Buscamos la fila que tiene el mayor valor de 
  # informacion mutua. Se obtiene como raiz la variable que condiciona en dicho valor
  root = ifelse(root==0,nrow(MI),root)# Al hacer el modulo, 0=root=nrow
  return(root)
}

# Funcion para aprender un tan-------------------------------------
fit_tan <- function(target, data, root = NULL, distributions = NULL,
                    mutualInfoCond = NULL,fit.args=NULL,all = F, parallel = F){
  
  if(is.null(target)){
   stop(paste0("Variable ", taget, " must be indicated")) 
  }
  
  if(!(target%in%colnames(data))){
    stop(paste0("Variable ", target, " must be a name of a variable of data")) 
  }
  if(is.null(fit.args)){
    fit.args=fit.args.null(fit.args)
  }
  # browser()
  # Obtencion de la informacion mutua y distribuciones-------------------
  # Obtencion de la estructura--------------------------------------
  
  if(is.null(mutualInfoCond)){
    MI_dist = mutual_information_tan(data,target = target,fit.args = fit.args,
                                     root = root,parallel = parallel)
    mutualInfoCond = MI_dist$MI
    # distributions = MI_dist$distributions
    root = MI_dist$root
  }
  
  adj_mat = fit_tan_structure(mutualInfoCond,target,root)
  var_i = colnames(data)[1]
  distTarget = distVar(data = data,target,fit.args)
  varPred = setdiff(colnames(data),target)
  parents = sapply(varPred,function(var_i){
    names(which(adj_mat[,var_i]==1))
  })
  if(parallel){
    library(doParallel)
    cores = detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)
    # browser()
    distributions = foreach(var_i = varPred,.packages = "MoTBFs",
                            .inorder = F)%dopar%{
      source("mutualInformation.R")
      distX_cond(data,var_i,parents[[var_i]],fit.args)
    }
    stopCluster(cl)
  }else{
    # browser()
    distributions = lapply(varPred,function(var_i){
      distX_cond(data,var_i,parents[[var_i]],fit.args)
    })
  }
  distributions = append(distributions,list(distTarget))
  
  bn = MoTBFs:::getFormatedBN(distributions)
  # dag = getDAG(bn)
  # graphviz.plot(dag)
  if(all){# Devolver informacion mutua y distribicones tambien
    return(list(bn=bn,mutualInfoCond=mutualInfoCond,distributions=distributions))
  }else{
    return(bn)
  }
}
## Example 1----------------------------------------------------

# set.seed(14)
# C = as.factor(sample(c("c1","c2"),100,T,c(0.4,0.6)))
# tbC = as.vector(table(C))
# set.seed(1452)
# X = c()
# X[C =="c1"] = rbeta(tbC[1],2,3)
# set.seed(548)
# X[C =="c2"] = 2+rbeta(tbC[2],2,3)
# set.seed(841)
# Y = as.factor(sample(c("y1","y2"),100,T))
# set.seed(74)
# Z = rnorm(100)
# data = data.frame(X,Y,Z,C)
# 
# target = "C"
# mutualInfoCond = NULL
# root=NULL
# distributions=NULL
# fit.args=NULL

# fit_tan(target = target,data = data,root=root)

## Example 2-----------------------------------------------------
# set.seed(14)
# C = as.factor(sample(c("c1","c2"),100,T,c(0.4,0.6)))
# tbC = as.vector(table(C))
# set.seed(1452)
# X = c()
# X[C =="c1"] = rbeta(tbC[1],2,3)
# set.seed(548)
# X[C =="c2"] = 2+rbeta(tbC[2],2,3)
# set.seed(841)
# Y = as.factor(sample(c("y1","y2"),100,T))
# set.seed(74)
# Z = rnorm(100)
# A= c()
# set.seed(12)
# A[C =="c1"] = sample(c("a1","a2"),tbC[1],T,c(0.7,0.3))
# set.seed(174)
# A[C =="c2"] = sample(c("a1","a2"),tbC[2],T,c(0.25,0.75))
# A = as.factor(A)
# target = "C"
# mutualInfoCond = NULL
# root=NULL
# distributions=NULL
# fit.args=fit.args.null(NULL)
# data = data.frame(X,Y,Z,C)
# data_new = data.frame(X,Y,Z,A,C)
# order = 1:4
# varsNew = "A"
# bn_old = fit_tan(target,data = data,all = T)
# mutualInfoCond_0 = bn_old$mutualInfoCond
# distributions_0 = bn_old$distributions
# # function(target, data_new, varsNew, root = NULL, distributions_0 = NULL,
# #          mutualInfoCond_0,fit.args=NULL)
# data_new
# MI_dist = update_mutual_information_tan(data_new = data_new,varsNew = varsNew,
#                                         target = target,
#                                         distributions0 = distributions_0,
#                                         mutualInfoCond_0 = mutualInfoCond_0,
#                                         fit.args = fit.args)
# mutualInfoCond = MI_dist$MI
# distributions = MI_dist$distributions
# root = MI_dist$root
# bnAll = fit_tan(target = target,data = data_new,root = root,
#                 distributions = distributions,
#                 mutualInfoCond = mutualInfoCond,
#                 all=T)
# dag=getDAG(bnAll$bn)
# graphviz.plot(dag)