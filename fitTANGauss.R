# Funciones para crear el TAN

# Función para el algoritmo de chow-liu--------------------
## Función para el arbol maximal--------------------
prim_maximal_g <- function(adj_mat,root = 1,pred_disc = character()) {
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
    if(all(selected_nodes[pred_disc])){
      option_nodes = names(which(!selected_nodes))
    }else{
      option_nodes = names(which(!selected_nodes[pred_disc]))
    }
    # Encontrar la arista de mayor peso que conecta un nodo seleccionado con uno no seleccionado
    for (j in names(which(selected_nodes))) {
      for (k in option_nodes) {
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
fit_tan_structure = function(mutualInfoCond,target,root=1,pred_disc){
  # browser()
  # Matriz de adjacencia
  adj_mat <- prim_maximal_g(mutualInfoCond,root = root)
  # nombre de las variables
  varTAN <- c(row.names(adj_mat),target)
  # anadimos la clase
  adj_mat <- rbind(adj_mat,1)# Aristas del target a todos los predictores
  adj_mat <- cbind(adj_mat,0)
  rownames(adj_mat)<- varTAN
  colnames(adj_mat)<- varTAN
  
  return(adj_mat)
}

# funcion para calcular la raiz de un TAN---------------------
fit_root = function(MI,pred_disc = character()){
  if(length(pred_disc)!=0){
    MI = MI[pred_disc,pred_disc,drop = F]
  }
  root = which.max(MI)%%nrow(MI)
  root = ifelse(root==0,nrow(MI),root)# Al hacer el modulo, 0=root=nrow
  root = colnames(MI)[root]
  
  return(root)
}

# Funcion para aprender un tan discreto CG--------------------------------
fit_tan_g <- function(target, data, root = NULL, distributions = list(),
                      mutualInfoCond = NULL, all = F){
  
  if(is.null(target)){
    stop(paste0("Variable ", taget, " must be indicated"))
  }
  
  if(!(target%in%colnames(data))){
    stop(paste0("Variable ", target, " must be a name of a variable of data"))
  }
  # browser()
  # Obtencion de la informacion mutua y distribuciones-------------------
  # Obtencion de la estructura--------------------------------------
  # hay variables discretas
  varPred = setdiff(colnames(data),target)
  pred_cont = sapply(varPred,function(variable){
    is.numeric(data[1,variable])
  })
  pred_cont = names(pred_cont)[which(pred_cont)]
  pred_disc = setdiff(varPred,pred_cont)
  # Obtenemos si fuese necesario la raiz
  if(is.null(root)||!(root%in%pred_disc)){
    root = fit_root(mutualInfoCond,pred_disc)
  }
  
  # Obtenemos si fuese necesario la raiz
  if(is.null(root)||!(root%in%pred_disc)){
    root = fit_root(mutualInfoCond,pred_disc)
  }
  # browser()
  adj_mat = fit_tan_structure(mutualInfoCond,target,root)
  dag = empty.graph(colnames(data))
  amat(dag)<-adj_mat
  bn = bn.fit(dag,data = data)
  if(all){# Devolver informacion mutua y distribicones tambien
    return(list(bn=bn,mutualInfoCond=mutualInfoCond,distributions=distributions))
  }else{
    return(bn)
  }
}

## Condicional Gaussiano--------------------------------------------
# X = varX
# data = datos[,c(target,X)]
mi_cond_gauss = function(X,target,data){
  levels_target = levels(data[[target]])
  varX = var(data[[X]])
  varX_t = sapply(levels_target, function(i){
    var(data[data[[target]]==i,X])
  })
  prob_target = as.vector(prop.table(table(data[[target]])))
  return(1/2*(log(varX)-sum(log(varX_t)*prob_target)))
}
pVar_C = function(discVar_C,pC){
  var = discVar_C$node
  pVar_C = data.frame(C = pC[["C"]],
                      mu = as.vector(discVar_C$coefficients[1,]),
                      sd = as.vector(discVar_C$sd))

  return(pVar_C)
}

cond_mi_cont = function(X,Y,target,data){
  levels_target = levels(data[[target]])
  prob_target = as.vector(prop.table(table(data[[target]])))
  
  corXY_t = sapply(levels_target, function(i){
    cor(data[data[[target]]==i,X],data[data[[target]]==i,Y])
  })
  return(-1/2*sum(prob_target*log((1-corXY_t^2))))
}

cond_mi_cont_disc = function(X,Y,target,data){
  browser()
  levels_target = levels(data[[target]])
  prob_target = as.vector(prop.table(table(data[[target]])))
  varX_t = sapply(levels_target, function(i){
    var(data[data[[target]]==i,X])
  })
  levels_Y = levels(data[[Y]])
  prob_joint = prop.table(table(data[,c(target,Y)]))
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
  # Posibles errores de redondeo o varianzas de 0
  var_Y_T = ifelse(var_Y_T<10^-6,10^-6,var_Y_T)
  return(1/2*(sum(prob_target*log(varX_t))-sum(prob_joint*log(var_Y_T),na.rm = T)))
}


MI_tan_gauss = function(data,target){
  
  # Información mutua entre discretas
  varPred = setdiff(colnames(data),target)
  MI = matrix(0,ncol = length(data)-1,nrow = length(data)-1,
              dimnames = list(varPred,varPred))
  for(Y in varPred){
    ff = function(X,Y,target,data){
      if(is.numeric(data[[X]])){
        if(is.numeric(data[[Y]])){
          return(cond_mi_cont(X,Y,target,data))
        }else{
          return(-1)
        }
      }else{
        if(is.numeric(data[[Y]])){
          return(cond_mi_cont_disc(Y,X,target,data))
        }else{
          return(infotheo::condinformation(data[[X]],data[[Y]],data[[target]]))
        }
      }
    }
    Mis = sapply(X = varPred, function(i){
      if(i==Y){
        return(0)
      }else{
        return(ff(i,Y,target,data[,c(i,Y,target)]))
      }
    })
    
    MI[names(Mis),Y] = Mis
  }
  return(MI)
}
