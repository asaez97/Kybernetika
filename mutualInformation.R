# Funciones de control------------------------------------------
fit.args.null = function(fit.args){
  if(is.null(fit.args$numIntervals)){
    fit.args$numIntervals=4
  }
  if(is.null(fit.args$POTENTIAL_TYPE)){
    fit.args$POTENTIAL_TYPE="MOP"
  }
  if(is.null(fit.args$scale)){
    fit.args$scale=F
  }
  if(is.null(fit.args$maxParam))
    fit.args$maxParam = 7
  return(fit.args)
}
# Funciones para calcular distribuciones---------------------------------------
distX_cond = function(data,varX, varsCond, fit.args =NULL){
  fit.args = fit.args.null(fit.args)
  # browser()
  varType = ifelse(is.numeric(data[,varX]),"Continuous","Discrete")
  distX_C = do.call(conditionalMethod,
                    append(fit.args,
                           list(data=data,nameParents=varsCond,
                                nameChild=varX)))
  distX_C = list(Child=varX,functions=distX_C,varType = varType)
  
  return(distX_C)
}

distVar = function(data,varName,fit.args){
  fit.args = fit.args.null(fit.args)
  if(is.numeric(data[,varName])){
    fit.args$x=data[,varName]
    fit.args$numIntervals=NULL
    distVar = do.call(univMoTBF,fit.args)
    distVar = list(Child=varName,functions=list(Px = distVar),# Para que este en el formato correcto
                  varType = "Continuous")
  }else{
    distVar = probDiscreteVariable(data[,varName])
    # Almacenamos la distribucion para devolverla
    distVar = list(Child=varName,functions=list(distVar),varType = "Discrete")
  }
  return(distVar)
}
# Funciones para integrar MoTBFs-----------------------------------------------
# Funcion para integrar en X
integrateX = function(CPD_X_,CPD_X,X){
  # browser()
  integrateX = c()
  for(i in 1:length(X)){
    domain = unlist(X[i])
    # tranformamos en funcion
    CPD_X_i = as.function(CPD_X_[i][[1]])
    CPD_Xi = as.function(CPD_X[i][[1]])
    # Calculo del integrando
    # integrando = function(x){
    #   xx = seq(domain[1],domain[2],10^-4)
    #   minn = max(10^-5,min(CPD_X_i(xx),CPD_Xi(xx)))
    #   return(CPD_X_i(x)*(log(CPD_X_i(x)/minn)-log(CPD_Xi(x)/minn)))
    #   # return(CPD_X_i(x)*(log(CPD_X_i(x)/CPD_Xi(x))))
    # }
    # integrateX[i]=integrate(integrando,domain[1],domain[2])$value
    integrando1 = function(x){
      # xx = seq(domain[1],domain[2],10^-4)
      # minn = max(10^-5,min(CPD_X_i(xx),CPD_Xi(xx)))
      return(CPD_X_i(x)*log(CPD_X_i(x)))
      # return(CPD_X_i(x)*(log(CPD_X_i(x)/CPD_Xi(x))))
    }
    integrando2 = function(x){
      # xx = seq(domain[1],domain[2],10^-4)
      # minn = max(10^-5,min(CPD_X_i(xx),CPD_Xi(xx)))
      return(CPD_X_i(x)*log(CPD_Xi(x)))
      # return(CPD_X_i(x)*(log(CPD_X_i(x)/CPD_Xi(x))))
    }
    integrateX[i]=tryCatch(integrate(integrando1,domain[1],domain[2])$value-
      integrate(integrando2,domain[1],domain[2])$value,
      error = function(e) {
        cat("Error:", e$message, "\n")
        0},warning = function(w) {
          cat("Advertencia:", w$message, "\n")
          0})
  }
  return(integrateX)
}

integrateVar = function(CPD,domain){
  
  integrateVar=c()
  i = 1
  for(i in 1:length(domain)){
    CPD_i = CPD[i][[1]]
    domain_i = unlist(domain[i])
    integrateVar[i] = integralMoTBF(CPD_i,domain_i[1],domain_i[2])
  }
  return(integrateVar)
}


# Informacion mutua para dos variables discretas-----------------

# fit.args =  Argumentos para la calcular las distribuciones MoTBFs
# Vease argumentos de ConditionalMethod. Si es null, emplea los
# argumentos por defecto de esa funcion con ademas numIntervals=4,
# POTENTIAL_TYPE = "MOP",scale = F
# Se puede introduccir la distribucion de C obtenida empleando la
# la funcion probDiscreteVariable()
# se puede introduccir tambien las distribuciones de Y|C, de X|C y de X| Y,C,
# con X la variable continua e Y la variable discreta. Todas aprendidas con 
# conditionalMethod
mut_information_dis = function(data, nameVars=NULL, fit.args =NULL,
                               distC = NULL,distX_C=NULL,distX = NULL){
  # Chequeamos los argumentos para fijar las Mops------------
  fit.args = fit.args.null(fit.args)
  # Control nameVars
  if(!is.null(nameVars)){
    data = data[,nameVars]
  }
  if(length(data)!=2){
    stop("data must have 2 variables")
  }else{
    nameVars = colnames(data)
  }
  
  numLevels = sapply(data,nlevels)
  if(any(numLevels<=1)){
    stop("Some variable is not discrete with more than one value")
  }
  # browser()
  X = nameVars[1]
  C = nameVars[2]
  distributions = list()
  nDist=0
  ## Distribucion de X|C------------------------------------
  if(is.null(distX_C)){
    # distX_C = conditionalMethod(data = data, nameParents = c("C"),
    #                                 nameChild = "X",numIntervals = 4,
    #                                 POTENTIAL_TYPE = "MOP",scale = scale)
    distX_C = distX_cond(data = data,varX = X,
                         varsCond = C,fit.args = fit.args)
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    distributions[[nDist]] = distX_C
    names(distributions)[nDist]=paste0(X,"|",C)
  }
  ## Nos quedamos solo con la CPD
  distX_C = MoTBFs:::getFormatedBN(list(distX_C))[[1]]$functions
  colnames(distX_C)[3]=paste0("CPD_",X,"|",C)
  
  # Distribucion de C----------------------------
  if(is.null(distC)){
    distC = distVar(data,varName = C,fit.args = fit.args)
    # Almacenamos la distribucion para devolverla
    nDist=nDist+1
    distributions[[nDist]] = distC
    names(distributions)[nDist]=C
  }
  ## Nos quedamos solo con la CPD
  distC =  MoTBFs:::getFormatedBN(list(distC))[[1]]$functions
  
  # Distribucion de X----------------------------
  if(is.null(distX)){
    distX = distVar(data,varName = X,fit.args = fit.args)
    # Almacenamos la distribucion para devolverla
    nDist=nDist+1
    distributions[[nDist]] = distX
    names(distributions)[nDist]=X
  }
  ## Nos quedamos solo con la CPD
  distX =  MoTBFs:::getFormatedBN(list(distX))[[1]]$functions
  
  # Arbol para los calculos
  arbol = merge.data.frame(distX_C,distX)
  arbol = merge.data.frame(arbol,distC)
  
  # Calculo de la informacion para cada rama
  arbol$sumX = arbol[,paste0("CPD_",X,"|",C)]*
    (log(arbol[,paste0("CPD_",X,"|",C)])-log(arbol[,paste0("CPD_",X)]))
  # Compute I(X,Y)
  MI = sum(arbol$sumX*arbol[,paste0("CPD_",C)])
  return(MI)
}
## Examples-----------
# varNames = NULL
# set.seed(1252)
# X = as.factor(sample(c("x1","x2","x3"),100,T,prob = c(0.3,0.4,0.3)))
# set.seed(1541)
# Y = as.factor(sample(c("y1","y2"),100,T,prob = c(0.4,0.6)))
# data = data.frame(X,Y)
# fit.args=NULL
# distX_C=NULL
# mut_information_dis(data = data)
# infotheo::multiinformation(data)

# Informacion mutua para variable X continua e Y discreta-----------------------
# varNames = c("mcg","lip")
# data = data2[,varNames]
mut_information_cont_dis = function(data,nameVars=NULL, fit.args =NULL,
                                    distC = NULL,distX_C=NULL,distX = NULL){
  fit.args = fit.args.null(fit.args)
  # browser()
  if(!is.null(nameVars)){
    data = data[,nameVars]
    if(is.null(data)){
      stop(paste0("data must have",nameVars,"as colnames",collapse = " "))
    }
  }else if(length(data)!=2){
    stop("data must have 2 variables")
  }else{
    nameVars = colnames(data)
  }
  numLevels = sapply(data,nlevels)
  if(all(numLevels)!=0){
    stop("There is not a continuous variable")
  }
  if(sum(numLevels)<=1){
    stop("There is not a discrete variable")
  }
  # browser()
  # Calculamos las distribuciones
  X = colnames(data)[numLevels==0]# Continua
  C = colnames(data)[numLevels!=0]# Discreta
  distributions = list()
  nDist=0
  # Distribucion de X----------------------------
  if(is.null(distX)){
    distX = distVar(data,varName = X,fit.args = fit.args)
    # Almacenamos la distribucion para devolverla
    nDist=nDist+1
    distributions[[nDist]] = distX
    names(distributions)[nDist]=X
  }
  ## Nos quedamos solo con la CPD
  distX =  MoTBFs:::getFormatedBN(list(distX))[[1]]$functions
  
  # Distribucion de C----------------------------
  if(is.null(distC)){
    distC = distVar(data,varName = C,fit.args = fit.args)
    # Almacenamos la distribucion para devolverla
    nDist=nDist+1
    distributions[[nDist]] = distC
    names(distributions)[nDist]=C
  }
  ## Nos quedamos solo con la CPD
  distC =  MoTBFs:::getFormatedBN(list(distC))[[1]]$functions
  ## Distribucion de X|C------------------------------------
  if(is.null(distX_C)){
    # distX_C = conditionalMethod(data = data, nameParents = c("C"),
    #                                 nameChild = "X",numIntervals = 4,
    #                                 POTENTIAL_TYPE = "MOP",scale = scale)
    distX_C = distX_cond(data = data,varX = X,
                         varsCond = C,fit.args = fit.args)
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    distributions[[nDist]] = distX_C
    names(distributions)[nDist]=paste0(X,"|",C)
  }
  ## Nos quedamos solo con la CPD
  distX_C = MoTBFs:::getFormatedBN(list(distX_C))[[1]]$functions
  colnames(distX_C)[3]=paste0("CPD_",X,"|",C)
  
  
  # Calculo de la Informacion Mutua-------------------
  
  # Construimos el arbol para los calculos
  arbol = cbind(distX_C,distX[,2,drop=F])
  # Incluimos la probabilidad de C
  arbol = merge.data.frame(arbol,distC)
  # Calculamos las diferentes integrales en X
  arbol$integrateX = integrateX(arbol[,paste0("CPD_",X,"|",C)],
                                arbol[,paste0("CPD_",X)],
                                arbol[,X])
  MI = sum(arbol[,paste0("CPD_",C)]*arbol$integrateX)
  return(MI)
}


# ## Examples-----------
# varNames = NULL
# set.seed(1252)
# X = rbeta(100,2,4)
# set.seed(1541)
# Y = as.factor(sample(c("y1","y2"),100,T,prob = c(0.4,0.6)))
# data = data.frame(X,Y)
# fit.args=NULL
# distX_C=NULL
# distC=NULL
# distX=NULL
# mut_information_cont_dis(data)

mut_information_dis_cont = function(data,nameVars=NULL, fit.args =NULL,
                                    distC = NULL,distC_X=NULL,distX = NULL){
  fit.args = fit.args.null(fit.args)
  # browser()
  if(!is.null(nameVars)){
    data = data[,nameVars]
    if(is.null(data)){
      stop(paste0("data must have",nameVars,"as colnames",collapse = " "))
    }
  }else if(length(data)!=2){
    stop("data must have 2 variables")
  }else{
    nameVars = colnames(data)
  }
  numLevels = sapply(data,nlevels)
  if(all(numLevels)!=0){
    stop("There is not a continuous variable")
  }
  if(sum(numLevels)<=1){
    stop("There is not a discrete variable")
  }
  # browser()
  # Calculamos las distribuciones
  X = colnames(data)[numLevels==0]# Continua
  C = colnames(data)[numLevels!=0]# Discreta
  distributions = list()
  nDist=0
  # Distribucion de X----------------------------
  if(is.null(distX)){
    distX = distVar(data,varName = X,fit.args = fit.args)
    # Almacenamos la distribucion para devolverla
    nDist=nDist+1
    distributions[[nDist]] = distX
    names(distributions)[nDist]=X
  }
  ## Nos quedamos solo con la CPD
  distX =  MoTBFs:::getFormatedBN(list(distX))[[1]]$functions
  
  # Distribucion de C----------------------------
  if(is.null(distC)){
    distC = distVar(data,varName = C,fit.args = fit.args)
    # Almacenamos la distribucion para devolverla
    nDist=nDist+1
    distributions[[nDist]] = distC
    names(distributions)[nDist]=C
  }
  ## Nos quedamos solo con la CPD
  distC =  MoTBFs:::getFormatedBN(list(distC))[[1]]$functions
  ## Distribucion de X|C------------------------------------
  if(is.null(distC_X)){
    # distX_C = conditionalMethod(data = data, nameParents = c("C"),
    #                                 nameChild = "X",numIntervals = 4,
    #                                 POTENTIAL_TYPE = "MOP",scale = scale)
    distC_X = distX_cond(data = data,varX = C,
                         varsCond = X,fit.args = fit.args)
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    distributions[[nDist]] = distC_X
    names(distributions)[nDist]=paste0(C,"|",X)
  }
  ## Nos quedamos solo con la CPD
  distC_X = MoTBFs:::getFormatedBN(list(distC_X))[[1]]$functions
  colnames(distC_X)[3]=paste0("CPD_",C,"|",X)
  
  
  # Calculo de la Informacion Mutua-------------------
  # browser()
  # Construimos el arbol para los calculos
  arbol = merge.data.frame(distC_X,distC)
  # Incluimos la probabilidad de C
  arbol = merge.data.frame(arbol,distX[,2,drop=F])
  # Calculamos las diferentes integrales en X
  arbol$sumC = log(arbol[[paste0("CPD_",C,"|",X)]]/arbol[[paste0("CPD_",C)]])*
    arbol[[paste0("CPD_",C,"|",X)]]
  arbol$integrateX = integrateVar(arbol[[paste0("CPD_",X)]],domain = arbol[[X]])
  MI = sum(arbol$sumC*arbol$integrateX)
  return(MI)
}


# Información mutua entre dos variables continuas----------------------------
# data = data[,1:2]
mut_information_cont = function(data,varNames=NULL){
  # Comprobaciones
  
  if(!is.null(varNames)){
    data = data[,varNames]
    if(is.null(data)){
      stop(paste0("data must have",varNames,"as colnames",collapse = " "))
    }
  }else if(length(data)!=2){
    stop("data must have 2 variables")
  }else{
    varNames = colnames(data)
  }
  numLevels = sapply(data,nlevels)
  if(sum(numLevels)!=0){
    stop("There is a non numerical variable")
  }
  
  # browser()
  # Calculamos las distribuciones
  X = data[,1]
  Y = data[,2]
  
  distX = univMoTBF(X,POTENTIAL_TYPE = "MOP",scale = T)
  px = as.function(distX)
  distY = univMoTBF(Y,POTENTIAL_TYPE = "MOP",scale = T)
  py = as.function(distY)
  
  distXgivenY =  conditionalMethod(data = data, nameParents = varNames[2],
                                   nameChild = varNames[1],numIntervals = 4,
                                   POTENTIAL_TYPE = "MOP",scale = T)
  # hist(X,probability = T)
  # plot(distX,add=T)
  # hist(X[Y=="b"],probability = T)
  # plot(distXgivenY[[2]]$Px,add=T)
  # Calculamos la información mutua
  
  # Bucle que pasa por cada particion
  result = 0
  i = 1
  for(i in 1:length(distXgivenY)){
    # Seleccionamos la distribucion X|Y \in y_i
    px_yi = as.function(distXgivenY[[i]]$Px)
    # DOminio de la variable X|Y\in I_i
    domain = sort(distXgivenY[[i]]$Px$Domain)
    # Intervalo de la particion de la variable Y
    intervals = distXgivenY[[i]]$interval
    
    # Calculo de la primera integral (X|Y \in y_i)
    integrando = function(x){
      return(px_yi(x)*log(px_yi(x)/px(x)))
    }
    integral = integrate(integrando,domain[1],domain[2])$value
    # Calculo integral Y*Integral1 si Y \in y_i
    result = result+integrate(py,intervals[1],intervals[2])$value*integral
  }
  return(result)
}

# Informacion mutua entre variables dada otra variable------


## Continua y discreta dada una discreta------------------------------

# fit.args =  Argumentos para la calcular las distribuciones MoTBFs
            # Vease argumentos de ConditionalMethod. Si es null, emplea los
            # argumentos por defecto de esa funcion con ademas numIntervals=4,
            # POTENTIAL_TYPE = "MOP",scale = F
# Se puede introduccir la distribucion de C obtenida empleando la
# la funcion probDiscreteVariable()
# se puede introduccir tambien las distribuciones de Y|C, de X|C y de X| Y,C,
# con X la variable continua e Y la variable discreta. Todas aprendidas con 
# conditionalMethod
cond_mut_information_cont_dis = function(data,className,corLap=1,
                                         varNames=NULL, fit.args =NULL,
                                         distC = NULL,distY_C=NULL,
                                         distX_YC=NULL,distX_C=NULL){
  if(is.null(fit.args$numIntervals)){
    fit.args$numIntervals=4
  }
  if(is.null(fit.args$POTENTIAL_TYPE)){
    fit.args$POTENTIAL_TYPE="MOP"
  }
  if(is.null(fit.args$scale)){
    fit.args$scale=F
  }
  # Comprobaciones
  # browser()
  if(!(className%in%colnames(data))){
    stop(paste0("data must have",className,"as colnames",collapse = " "))
  }
  # C = data[,className]
  if(!is.null(varNames)){
    data = data[,varNames]
    if(is.null(data)){
      stop(paste0("data must have",varNames,"as colnames",collapse = " "))
    }
  }
  if(length(data)!=3){
    stop("data must have 3 variables")
  }else{
    varNames = colnames(data)
  }
  # Eliminamos la variable clase del data.frame
  data2 = data[,colnames(data)!=className]
  numLevels = sapply(data2,nlevels)
  if(all(numLevels)!=0){
    stop("There is not a continuous variable")
  }
  if(sum(numLevels)<=1){
    stop("There is not a discrete variable")
  }
  
  varCont = colnames(data2)[numLevels==0]# X
  varDisc = colnames(data2)[numLevels!=0]# Y
  rm(data2)
  
  distributions = list()
  # Cuenta el numero nuevo de distribuciones aprendidas
  nDist = length(distributions)
  # browser()
  # # Calculamos las distribuciones
  # X = data2[,numLevels==0]# Continua
  # Y = data2[,numLevels!=0]# Discreta
  # data = data.frame(X,Y,C)
  
  
  # Distribucion de X|C,Y
  # distX_YC = conditionalMethod(data,c("C","Y"),"X",4,"MOP",scale=F)
  if(is.null(distX_YC)){
    distX_YC = do.call(conditionalMethod,
                       append(fit.args,
                              list(data=data,nameParents=c(className,varDisc),
                                   nameChild=varCont)))
    distX_YC = list(Child=varCont,functions=distX_YC,varType = "Continuous")
    nDist=nDist+1
    distributions[[nDist]] = distX_YC
    parents = sort(c(className,varDisc))
    names(distributions)[nDist]=paste0(varCont,"|",paste0(parents,collapse = ":"))
    
  }
  # Almacenamos la distribucion para devolverla
  
  ## Nos quedamos solo con la CPD
  distX_YC = MoTBFs:::getFormatedBN(list(distX_YC))[[1]]$functions
  colnames(distX_YC)[4]=paste0("CPD_",varCont,"|",varDisc,":",className)
  
  
  ## Distribucion de X|C------------------------------------
  if(is.null(distX_C)){
    # distX_C = conditionalMethod(data = data, nameParents = c("C"),
    #                                 nameChild = "X",numIntervals = 4,
    #                                 POTENTIAL_TYPE = "MOP",scale = scale)
    distX_C = do.call(conditionalMethod,
                      append(fit.args,
                             list(data=data,nameParents=className,
                                  nameChild=varCont)))
    distX_C = list(Child=varCont,functions=distX_C,varType = "Continuous")
    
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    distributions[[nDist]] = distX_C
    names(distributions)[nDist]=paste0(varCont,"|",className)
  }
  ## Nos quedamos solo con la CPD
  distX_C = MoTBFs:::getFormatedBN(list(distX_C))[[1]]$functions
  colnames(distX_C)[3]=paste0("CPD_",varCont,"|",className)
  
  # Distribucion de C----------------------------
  if(is.null(distC)){
    # distC = prop.table((table(data$C)+1))
    distC = probDiscreteVariable(data[,className])
    
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    # Almacenamos la distribucion para devolverla
    distC = list(Child=className,functions=list(distC),varType = "Discrete")
    distributions[[nDist]] = distC
    names(distributions)[nDist]=paste0(className)
  }
  ## Nos quedamos solo con la CPD
  distC =  MoTBFs:::getFormatedBN(list(distC))[[1]]$functions
  
  
  # Distribucion de Y|C--------------------------
  if(is.null(distY_C)){
    # Distribucion de Y dada C
    # distY_C = conditionalMethod(data,"C","Y",4,"MOP",scale = scale)
    distY_C = do.call(conditionalMethod,
                      append(fit.args,
                             list(data=data,nameParents=className,
                                  nameChild=varDisc)))
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    distY_C = list(Child=varDisc,functions=distY_C,varType = "Discrete")
    distributions[[nDist]] = distY_C
    names(distributions)[nDist]=paste0(varDisc,"|",className)
  }
  ## Nos quedamos solo con la CPD
  distY_C = MoTBFs:::getFormatedBN(list(distY_C))[[1]]$functions
  colnames(distY_C)[3]=paste0("CPD_",varDisc,"|",className)
  # Calculo de la Informacion Mutua-------------------
  
  # Construimos el arbol para los calculos
  arbol = merge.data.frame(distX_YC,distX_C)
  
  # Calculamos las diferentes integrales en X
  arbol$integrateX = integrateX(arbol[,paste0("CPD_",varCont,"|",varDisc,
                                              ":",className)],
                                arbol[,paste0("CPD_",varCont,"|",className)],
                                arbol[,varCont])
  
  # Incluimos la probabilidad de C
  arbol = merge.data.frame(arbol,distC)
  
  # Añadimos las probabilidades Y|C al arbol
  arbol = merge.data.frame(arbol,distY_C)
  

  MI = sum(arbol$integrateX*
             arbol[,paste0("CPD_",varDisc,"|",className)]*
             arbol[,paste0("CPD_",className)])
  
  return(list(MI=MI,distributions=distributions))
}
# 
# ### Examples-------
# data(ecoli,package = "MoTBFs")
# className = "lip"
# varNames = c("alm2","chg",className)
# scale = F
# data = MoTBFs:::check_data(ecoli[,-c(1,9)])
# corLap = 1
# fit.args = c()
# cond_mut_information_cont_dis(data,className,scale,varNames)



## Continuas dada una discreta------------------------------
# fit.args =  Argumentos para la calcular las distribuciones MoTBFs
# Vease argumentos de ConditionalMethod. Si es null, emplea los
# argumentos por defecto de esa funcion con ademas numIntervals=4,
# POTENTIAL_TYPE = "MOP",scale = F
# Se puede introduccir la distribucion de C obtenida empleando la
# la funcion probDiscreteVariable()
# se puede introducir tambien las distribuciones de Y|C, de X|C y de X| Y,C,
# con X e Y las variables continuas. Todas aprendidas con 
# conditionalMethod
cond_mut_information_cont = function(data,className,corLap=1,
                                     varNames=NULL, fit.args =NULL,
                                     distC = NULL,distY_C=NULL,
                                     distX_YC=NULL,distX_C=NULL){
  if(is.null(fit.args$numIntervals)){
    fit.args$numIntervals=4
  }
  if(is.null(fit.args$POTENTIAL_TYPE)){
    fit.args$POTENTIAL_TYPE="MOP"
  }
  if(is.null(fit.args$scale)){
    fit.args$scale=F
  }
  # Comprobaciones
  
  # browser()
  if(!(className%in%colnames(data))){
    stop(paste0("data must have",className,"as colnames",collapse = " "))
  }
  # C = data[,className]
  
  if(!is.null(varNames)){
    if(!all(varNames%in%colnames(data))){
      stop(paste0("data must have",varNames,"as colnames",collapse = " "))
    }
    data = data[,varNames]
  }
  
  if(length(data)!=3){
    stop("data must have 3 variables")
  }else{
    varNames = colnames(data)
  }
  
  # Eliminamos la variable clase del data.frame
  data2 = data[,colnames(data)!=className]
  
  numLevels = sapply(data2,nlevels)
  if(any(numLevels)!=0){
    stop("There is a discrete variable")
  }
  
  # browser()
  # Calculamos las distribuciones
  varCont1 = colnames(data2)[1]# X
  varCont2 = colnames(data2)[2]# Y
  rm(data2)
  
  distributions = list()
  # Cuenta el numero nuevo de distribuciones aprendidas
  nDist = length(distributions)
  
  # Distribucion de X|C,Y---------------------------------
  # distX_YC = conditionalMethod(data,c("C","Y"),"X",4,"MOP",scale=F)
  if(is.null(distX_YC)){
    distX_YC = do.call(conditionalMethod,
                       append(fit.args,
                              list(data=data,nameParents=c(className,varCont2),
                                   nameChild=varCont1)))
    # Almacenamos la distribucion para devolverla
    distX_YC = list(Child=varCont1,functions=distX_YC,varType = "Continuous")
    nDist=nDist+1
    distributions[[nDist]] = distX_YC
    parents = sort(c(className,varCont2))
    names(distributions)[nDist]=paste0(varCont1,"|",paste0(parents,collapse = ":"))
    
  }
  
  ## Nos quedamos solo con la CPD
  distX_YC = MoTBFs:::getFormatedBN(list(distX_YC))[[1]]$functions
  colnames(distX_YC)[4]=paste0("CPD_",varCont1,"|",varCont2,":",className)
  
  
  # Distribucion de X|C-----------------------------
  if(is.null(distX_C)){
    # distX_C = conditionalMethod(data = data, nameParents = c("C"),
    #                                 nameChild = "X",numIntervals = 4,
    #                                 POTENTIAL_TYPE = "MOP",scale = scale)
    distX_C = do.call(conditionalMethod,
                      append(fit.args,
                             list(data=data,nameParents=className,
                                  nameChild=varCont1)))
    distX_C = list(Child=varCont1,functions=distX_C,varType = "Continuous")
    
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    distributions[[nDist]] = distX_C
    names(distributions)[nDist]=paste0(varCont1,"|",className)
  }
  ## Nos quedamos solo con la CPD
  distX_C = MoTBFs:::getFormatedBN(list(distX_C))[[1]]$functions
  colnames(distX_C)[3]=paste0("CPD_",varCont1,"|",className)
  
  # Distribucion de C----------------------------
  if(is.null(distC)){
    # distC = prop.table((table(data$C)+1))
    distC = probDiscreteVariable(data[,className])
    
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    # Almacenamos la distribucion para devolverla
    distC = list(Child=className,functions=list(distC),varType = "Discrete")
    distributions[[nDist]] = distC
    names(distributions)[nDist]=paste0(className)
  }
  ## Nos quedamos solo con la CPD
  distC =  MoTBFs:::getFormatedBN(list(distC))[[1]]$functions
  
  
  # Distribucion de Y|C--------------------------
  if(is.null(distY_C)){
    # Distribucion de Y dada C
    # distY_C = conditionalMethod(data,"C","Y",4,"MOP",scale = scale)
    distY_C = do.call(conditionalMethod,
                      append(fit.args,
                             list(data=data,nameParents=className,
                                  nameChild=varCont2)))
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    distY_C = list(Child=varCont2,functions=distY_C,varType = "Continuous")
    distributions[[nDist]] = distY_C
    names(distributions)[nDist]=paste0(varCont2,"|",className)
  }
  ## Nos quedamos solo con la CPD
  distY_C = MoTBFs:::getFormatedBN(list(distY_C))[[1]]$functions
  colnames(distY_C)[3]=paste0("CPD_",varCont2,"|",className)
  
  # Calculo de la Informacion Mutua-------------------
  
  # Construimos el arbol para los calculos
  arbol = merge.data.frame(distX_YC,distX_C)
  
  # Calculamos las diferentes integrales en X
  arbol$integrateX = integrateX(arbol[,paste0("CPD_",varCont1,"|",varCont2,
                                              ":",className)],
                                arbol[,paste0("CPD_",varCont1,"|",className)],
                                arbol[,varCont1])
  
  # Incluimos la probabilidad de C
  arbol = merge.data.frame(arbol,distC)
  
  # Añadimos las probabilidades Y|C al arbol
  arbol = merge.data.frame(arbol,distY_C[,-2])# No incluir el dominio de Y
  
  arbol$integrateY = integrateVar(arbol[,paste0("CPD_",varCont2,"|",className)],
                                  arbol[,varCont2])
  
  # Calculamos la informacion mutua
  MI = sum(arbol$integrateX*arbol$integrateY*arbol[,paste0("CPD_",className)])
  return(list(MI=MI,distributions=distributions))
}


# ### Examples-------
# library(MoTBFs)
# data(ecoli,package = "MoTBFs")
# className = "lip"
# varNames = c("gvh","alm2",className)
# scale = F
# data = MoTBFs:::check_data(ecoli[,-c(1,9)])
# cond_mut_information_cont(data,className,scale,varNames)
# jointmotbf.fit(data[,varNames[1:2]])
# Pjoint = jointmotbf.fit(data[,varNames[1:2]])
# plot(Pjoint)


## Discretas dada una discreta------------------------------

###-----------------------------##
# Importante Condiciona la segunda variable del data.frame
###--------------------------###

# fit.args =  Argumentos para la calcular las distribuciones MoTBFs
# Vease argumentos de ConditionalMethod. Si es null, emplea los
# argumentos por defecto de esa funcion con ademas numIntervals=4,
# POTENTIAL_TYPE = "MOP",scale = F
# Se puede introduccir la distribucion de C obtenida empleando la
# la funcion probDiscreteVariable()
# se puede introducir tambien las distribuciones de Y|C, de X|C y de X| Y,C,
# con X e Y las variables discretas. Todas aprendidas con 
# conditionalMethod

cond_mut_information_dis = function(data,className,corLap=1,
                                    varNames=NULL, fit.args =NULL,
                                    distC = NULL,distY_C=NULL,
                                    distX_YC=NULL,distX_C=NULL){
  if(is.null(fit.args$numIntervals)){
    fit.args$numIntervals=4
  }
  if(is.null(fit.args$POTENTIAL_TYPE)){
    fit.args$POTENTIAL_TYPE="MOP"
  }
  if(is.null(fit.args$scale)){
    fit.args$scale=F
  }
  # Comprobaciones
  # browser()
  if(!(className%in%colnames(data))){
    stop(paste0("data must have",className,"as colnames",collapse = " "))
  }
  # C = data[,className]
  
  if(!is.null(varNames)){
    if(!all(varNames%in%colnames(data))){
      stop(paste0("data must have",varNames,"as colnames",collapse = " "))
    }
    data = data[,varNames]
  }
  
  if(length(data)!=3){
    stop("data must have 3 variables")
  }else{
    varNames = colnames(data)
  }
  
  # Eliminamos la variable clase del data.frame
  data2 = data[,colnames(data)!=className]
  
  numLevels = sapply(data2,nlevels)
  if(any(numLevels)==0){
    stop("There is a Discrete variable")
  }
  
  # browser()
  # Calculamos las distribuciones
  varDisc1 = colnames(data2)[1]# X
  varDisc2 = colnames(data2)[2]# Y
  rm(data2)
  
  distributions = list()
  # Cuenta el numero nuevo de distribuciones aprendidas
  nDist = length(distributions)
  
  
  # Distribucion de X|C,Y---------------------------------
  # distX_YC = conditionalMethod(data,c("C","Y"),"X",4,"MOP",scale=F)
  if(is.null(distX_YC)){
    distX_YC = do.call(conditionalMethod,
                       append(fit.args,
                              list(data=data,nameParents=c(className,varDisc2),
                                   nameChild=varDisc1)))
    # Almacenamos la distribucion para devolverla
    distX_YC = list(Child=varDisc1,functions=distX_YC,varType = "Discrete")
    nDist=nDist+1
    distributions[[nDist]] = distX_YC
    parents=paste0(sort(c(className,varDisc2)),collapse = ":")
    names(distributions)[nDist]=paste0(varDisc1,"|",parents)
    
  }
  
  ## Nos quedamos solo con la CPD
  distX_YC = MoTBFs:::getFormatedBN(list(distX_YC))[[1]]$functions
  colnames(distX_YC)[4]=paste0("CPD_",varDisc1,"|",varDisc2,":",className)
  
  
  # Distribucion de X|C-----------------------------
  if(is.null(distX_C)){
    # distX_C = conditionalMethod(data = data, nameParents = c("C"),
    #                                 nameChild = "X",numIntervals = 4,
    #                                 POTENTIAL_TYPE = "MOP",scale = scale)
    distX_C = do.call(conditionalMethod,
                      append(fit.args,
                             list(data=data,nameParents=className,
                                  nameChild=varDisc1)))
    distX_C = list(Child=varDisc1,functions=distX_C,varType = "Discrete")
    
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    distributions[[nDist]] = distX_C
    names(distributions)[nDist]=paste0(varDisc1,"|",className)
  }
  ## Nos quedamos solo con la CPD
  distX_C = MoTBFs:::getFormatedBN(list(distX_C))[[1]]$functions
  colnames(distX_C)[3]=paste0("CPD_",varDisc1,"|",className)
  
  # Distribucion de C----------------------------
  if(is.null(distC)){
    # distC = prop.table((table(data$C)+1))
    distC = probDiscreteVariable(data[,className])
    
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    # Almacenamos la distribucion para devolverla
    distC = list(Child=className,functions=list(distC),varType = "Discrete")
    distributions[[nDist]] = distC
    names(distributions)[nDist]=paste0(className)
  }
  ## Nos quedamos solo con la CPD
  distC =  MoTBFs:::getFormatedBN(list(distC))[[1]]$functions
  
  
  # Distribucion de Y|C--------------------------
  if(is.null(distY_C)){
    # Distribucion de Y dada C
    # distY_C = conditionalMethod(data,"C","Y",4,"MOP",scale = scale)
    distY_C = do.call(conditionalMethod,
                      append(fit.args,
                             list(data=data,nameParents=className,
                                  nameChild=varDisc2)))
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    distY_C = list(Child=varDisc2,functions=distY_C,varType = "Discrete")
    distributions[[nDist]] = distY_C
    names(distributions)[nDist]=paste0(varDisc2,"|",className)
  }
  ## Nos quedamos solo con la CPD
  distY_C = MoTBFs:::getFormatedBN(list(distY_C))[[1]]$functions
  colnames(distY_C)[3]=paste0("CPD_",varDisc2,"|",className)
  
  # Calculo de la Informacion Mutua-------------------
  
  # Construimos el arbol para los calculos
  arbol = merge.data.frame(distX_YC,distX_C)
  # Incluimos la probabilidad de C
  arbol = merge.data.frame(arbol,distC)
  
  # Añadimos las probabilidades Y|C al arbol
  arbol = merge.data.frame(arbol,distY_C)
  
  # Calculo del log(p(X|YC)/p(X|C))
  arbol$log = arbol[,paste0("CPD_",varDisc1,"|",varDisc2,":",className)]*
    (log(arbol[,paste0("CPD_",varDisc1,"|",varDisc2,":",className)])-
       log(arbol[,paste0("CPD_",varDisc1,"|",className)]))
  
  # Calculamos la informacion mutua
  MI = sum(arbol$log*
             arbol[,paste0("CPD_",varDisc2,"|",className)]*
             arbol[,paste0("CPD_",className)])
  return(list(MI=MI,distributions=distributions))
}


### Examples-------
# library(MoTBFs)
# set.seed(1)
# X = sample(letters[1:3],100,replace = T)
# 
# set.seed(21)
# Y = sample(c("1","2"),100,replace = T,prob=c(0.75,0.25))
# 
# set.seed(48)
# C = sample(c("C1","C2"),100,replace = T,prob=c(0.35,0.65))
# 
# data = data.frame(X,Y,C)
# className = "C"
# varNames = c("X","Y",className)
# scale = F
# data = MoTBFs:::check_data(data)
# 
# cond_mut_information_dis(data,className,scale,varNames)
# jointmotbf.fit(data[,varNames[1:2]])
# Pjoint = jointmotbf.fit(data[,varNames[1:2]])
# plot(Pjoint)

## Discreta y continua dada una discreta------------------------------
### Quien condiciona es la continua
# fit.args =  Argumentos para la calcular las distribuciones MoTBFs
# Vease argumentos de ConditionalMethod. Si es null, emplea los
# argumentos por defecto de esa funcion con ademas numIntervals=4,
# POTENTIAL_TYPE = "MOP",scale = F
# Se puede introduccir la distribucion de C obtenida empleando la
# la funcion probDiscreteVariable()
# se puede introduccir tambien las distribuciones de Y|C, de X|C y de X| Y,C,
# con X la variable continua e Y la variable discreta. Todas aprendidas con 
# conditionalMethod

cond_mut_information_dis_cont = function(data,className,corLap=1,
                                         varNames=NULL, fit.args =NULL,
                                         distC = NULL,distY_C=NULL,
                                         distY_XC=NULL,distX_C=NULL){
  if(is.null(fit.args$numIntervals)){
    fit.args$numIntervals=4
  }
  if(is.null(fit.args$POTENTIAL_TYPE)){
    fit.args$POTENTIAL_TYPE="MOP"
  }
  if(is.null(fit.args$scale)){
    fit.args$scale=F
  }
  # Comprobaciones
  # browser()
  if(!(className%in%colnames(data))){
    stop(paste0("data must have",className,"as colnames",collapse = " "))
  }
  
  if(!is.null(varNames)){
    data = data[,varNames]
    if(is.null(data)){
      stop(paste0("data must have",varNames,"as colnames",collapse = " "))
    }
  }
  if(length(data)!=3){
    stop("data must have 3 variables")
  }else{
    varNames = colnames(data)
  }
  # Eliminamos la variable clase del data.frame
  data2 = data[,colnames(data)!=className]
  numLevels = sapply(data2,nlevels)
  if(all(numLevels)!=0){
    stop("There is not a continuous variable")
  }
  if(sum(numLevels)<=1){
    stop("There is not a discrete variable")
  }
  
  varCont = colnames(data2)[numLevels==0]# X
  varDisc = colnames(data2)[numLevels!=0]# Y
  rm(data2)
  
  distributions = list()
  # Cuenta el numero nuevo de distribuciones aprendidas
  nDist = length(distributions)
  
  # browser()
  # Calculamos las distribuciones------------
  
  ## Distribucion de Y|C,X-----------------------
  if(is.null(distY_XC)){
    distY_XC = do.call(conditionalMethod,
                       append(fit.args,
                              list(data=data,nameParents=c(className,varCont),
                                   nameChild=varDisc)))
    distY_XC = list(Child=varDisc,functions=distY_XC,varType = "Discrete")
    nDist=nDist+1
    distributions[[nDist]] = distY_XC
    parents=paste0(sort(c(className,varCont)),collapse = ":")
    names(distributions)[nDist]=paste0(varDisc,"|",parents)
    
  }
  # Almacenamos la distribucion para devolverla
  
  ## Nos quedamos solo con la CPD
  distY_XC = MoTBFs:::getFormatedBN(list(distY_XC))[[1]]$functions
  colnames(distY_XC)[4]=paste0("CPD_",varDisc,"|",varCont,":",className)
  
  ## Distribucion de X|C------------------------------------
  if(is.null(distX_C)){
    # distX_C = conditionalMethod(data = data, nameParents = c("C"),
    #                                 nameChild = "X",numIntervals = 4,
    #                                 POTENTIAL_TYPE = "MOP",scale = scale)
    distX_C = do.call(conditionalMethod,
                      append(fit.args,
                             list(data=data,nameParents=className,
                                  nameChild=varCont)))
    distX_C = list(Child=varCont,functions=distX_C,varType = "Continuous")
    
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    distributions[[nDist]] = distX_C
    names(distributions)[nDist]=paste0(varCont,"|",className)
  }
  ## Nos quedamos solo con la CPD
  distX_C = MoTBFs:::getFormatedBN(list(distX_C))[[1]]$functions
  colnames(distX_C)[3]=paste0("CPD_",varCont,"|",className)
  
  
  # Distribucion de C----------------------------
  if(is.null(distC)){
    # distC = prop.table((table(data$C)+1))
    distC = probDiscreteVariable(data[,className])
    
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    # Almacenamos la distribucion para devolverla
    distC = list(Child=className,functions=list(distC),varType = "Discrete")
    distributions[[nDist]] = distC
    names(distributions)[nDist]=paste0(className)
  }
  ## Nos quedamos solo con la CPD
  distC =  MoTBFs:::getFormatedBN(list(distC))[[1]]$functions
  
  
  ## Distribucion de Y|C--------------------------
  if(is.null(distY_C)){
    # Distribucion de Y dada C
    # distY_C = conditionalMethod(data,"C","Y",4,"MOP",scale = scale)
    distY_C = do.call(conditionalMethod,
                      append(fit.args,
                             list(data=data,nameParents=className,
                                  nameChild=varDisc)))
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    distY_C = list(Child=varDisc,functions=distY_C,varType = "Discrete")
    distributions[[nDist]] = distY_C
    names(distributions)[nDist]=paste0(varDisc,"|",className)
  }
  ## Nos quedamos solo con la CPD
  distY_C = MoTBFs:::getFormatedBN(list(distY_C))[[1]]$functions
  colnames(distY_C)[3]=paste0("CPD_",varDisc,"|",className)
  
  
  # Construimos el arbol para los calculos
  arbol = merge.data.frame(distY_XC,distY_C)
  # Incluimos la probabilidad de C
  arbol = merge.data.frame(arbol,distC)
  #Incluimos la distribucion de X|C
  # no modificar el dominio de X
  arbol = merge.data.frame(arbol,distX_C[,-2])
  
  # p(y|xc)*(log(y|xc)-log(y|c))
  arbol$sumY = arbol[,paste0("CPD_",varDisc,"|",varCont,":",className)]*
    (log(arbol[,paste0("CPD_",varDisc,"|",varCont,":",className)])-
       log(arbol[,paste0("CPD_",varDisc,"|",className)]))
  
  # Calculamos las diferentes integrales en X
  arbol$integrateX = integrateVar(arbol[,paste0("CPD_",varCont,"|",className)],
                                  arbol[,varCont])
  # Calculo de la Informacion Mutua
  MI = sum(arbol$integrateX*arbol$sumY*arbol[,paste0("CPD_",className)])
  
  
  return(list(MI=MI,distributions=distributions))
}

### Example------------------
# library(MoTBFs)
# set.seed(1)
# X = rnorm(100)
# set.seed(21)
# Y = as.factor(sample(c("y1","y2","y3"),100,replace = T,prob=c(0.45,0.25,0.3)))
# 
# set.seed(48)
# C = as.factor(sample(c("C1","C2"),100,replace = T,prob=c(0.35,0.65)))
# 
# data = data.frame(X,Y,C)
# className = "C"
# varNames = c("X","Y",className)
# scale=F
# cond_mut_information_dis_cont(data,className,scale,varNames)

# Funcion para calcular la informacion mutua de una variable con el resto de un data.frame---------------------
mutual_info_df_Y = function(data,varCond,fit.args,distr=NULL){
  MI = c()
  varNames = setdiff(colnames(data),varCond)
  if(is.numeric(data[,varCond])){
    stop("No esta implementado MI con variable continuas como objetivo")
  }
  # var_i = varNames[3]
  for(var_i in varNames){
    if(is.numeric(data[,var_i])){
      MI_i = mut_information_cont_dis(data = data[,c(var_i,varCond)],
                                      fit.args = fit.args,distC = distr[[varCond]],
                                      distX_C = distr[[paste0(var_i,"|",varCond)]],
                                      distX = distributions[[var_i]])
      MI = c(MI,MI_i)
    }else{
      MI_i = mut_information_dis(data = data[,c(var_i,varCond)],
                                 distC = distr[[varCond]],
                                 distX_C = distr[[paste0(var_i,"|",varCond)]],
                                 distX = distributions[[var_i]])
      MI = c(MI,MI_i)
    }
  }
  names(MI)=varNames
  return(MI)
}

cond_mutual_info= function(X,Y,C,data,distC,distX_C,distY_C,fit.args){
  # browser()
  if(X==Y){
    return(list(MI=0,name = Y))
  }else{
    if(is.null(distC)){
      
    }
    if(is.null(distX_C)){
      
    }
    
    if(is.null(distX_C)){
      
    }
    if(distX_C$varType=="Continuous"){
      if(distY_C$varType=="Continuous"){
        result = cond_mut_information_cont(data = data,className = target,
                                           varNames = c(X,Y,target),
                                           fit.args = fit.args,
                                           distC = distC,distY_C = distY_C,
                                           distX_C = distX_C)
      }else{# Quien condiciona es discreto
        result = cond_mut_information_cont_dis(data = data,className = target,
                                               varNames = c(X,Y,target),
                                               fit.args = fit.args,
                                               distC = distC,distY_C = distY_C,
                                               distX_C = distX_C)
      }
    }else{
      # La distribucion de X es discreta
      if(distY_C$varType=="Continuous"){
        result = cond_mut_information_dis_cont(data = data,className = target,
                                               varNames = c(X,Y,target),
                                               fit.args = fit.args,
                                               distC = distC,distY_C = distX_C,
                                               distX_C = distY_C)
      }else{# Quien condiciona es discreto
        result = cond_mut_information_dis(data = data,className = target,
                                          varNames = c(X,Y,target),
                                          fit.args = fit.args,
                                          distC = distC,distY_C = distY_C,
                                          distX_C = distX_C)
      }
    }
  }
  return(MI = result$MI)
}
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
# 
# varCond="C"
# fit.args =NULL
# distr=NULL
# mutual_info_df_Y(data = data,varCond = varCond,fit.args = fit.args)
