############################################################
############ Fonctions de bases à utiliser dans R ####################
############ Ex: "VEM_last" ####################
############ Inférence du modele par VEM ####################
############################################################

##########################################################################
# Inference du model SBM avec melange par VEM
VEM <- function(S, K, niter=100, epsilon_tau=1e-4, epsilon_eta=1e-4, verbose=FALSE, explorFact=1.5){
   # niter=100; epsilon_tau=1e-4; epsilon_eta = 1e-4; verbose = FALSE
  
##########################################################################
#gere les Scores d'entrées matrice ou vecteur et les rend en vecteur
  
  if(is.null(dim(S)))
     {
     N <- length(S)
     n <-(1+sqrt(1+8*N))/2
     transfo_indices <- indices(n)
     mat_S <- vect_mat_low(S)
     vec_S <- S
     }else{
        n=length(S[1,])
        N=(n-1) * n / 2
        transfo_indices <- indices(n)
        mat_S <- S
        vec_S <- mat_vect_low(S)
     }
   
   A=array(0,dim=c(N,K,K))
   
##############################################################################  
  
  #initialisation des paramètres eta et G (G a deux classes arete ou pas)
  param_gm <- Mclust(vec_S, G=2, verbose=FALSE) 
  eta_init<- param_gm$z
  g_init <- param_gm$classification-1
  
  #par la moyenne reclassement des scores avec et sans aretes
  mean_by_class <- vapply(0:1,function(g){mean(vec_S[g_init==g])},1)
  if(mean_by_class[2]  < mean_by_class[1] ) { 
    eta_init <- eta_init[,c(2,1)]
    g_init <- 1 - g_init
  }

  eta0 <- array(rep(eta_init[,1],K*K),c(N,K,K)) #ok
  eta1 <- array(rep(eta_init[,2],K*K),c(N,K,K)) 
  
  #initialisation des tau et  Pi comme la moyenne des tau de chaque classe
  param_sbm <- BM_bernoulli(membership_type="SBM_sym", adj=vect_mat_low( g_init), 
                            plotting='', verbosity=0, ncores=1, exploration_factor=explorFact)
  param_sbm$estimate()
  tau_init <- param_sbm$memberships[[K]]$Z
  tau_hat <- tau_init
  Pi_hat = colMeans(tau_hat)

  #init des normales, distribution des normales
  phi0 <-phi1<-  rep(0,2)

  #init du vecteur de borne inf pour plot
  vec_BI <- rep(0, 3*niter)

  #itération 
  diff = 2*epsilon_tau
  t<-0
  
  while ((t < niter) & (diff>epsilon_tau)){cat(t, '')
    t <- t+1
    #Etape M:
    
    #Calcul gamma:
    Gamma_hat_num <- matrix(nrow = K, ncol = K)
    Gamma_hat_num <- sapply(1:K, function(k){ sapply(1:K, function(l){ 
      t(tau_hat[,k]) %*% vect_mat_low(eta1[,k,l]) %*% tau_hat[,l] 
      })})
    Gamma_hat_denum <- (t(tau_hat) %*% (matrix(1,n,n) - diag(1,n)) %*% tau_hat)
    Gamma_hat <- Gamma_hat_num / Gamma_hat_denum
    
    #calcul de Pi:
    Pi_hat <- colMeans(tau_hat)
    
    #calcul de la borne_inf pour control:
    vec_BI[(3*t)-2] =  borne_inf(tau_hat,Pi_hat,A,eta0,eta1)

    #calcul de phi les parametres des densités
    Psi0_hat = Psi1_hat = rep(0, N)    
    invisible(sapply(1:(n-1), function(j){sapply((j+1):n, function(i){
      ij = which((transfo_indices[, 1]==i) & (transfo_indices[, 2]==j))
      # cat(i, j, ij, transfo_indices[ij, ], '\n')
      Psi0_hat[ij] <<-  (t(tau_hat[i,]) %*% eta0[ij, , ] %*% tau_hat[j,])[1, 1]
      Psi1_hat[ij] <<-  (t(tau_hat[i,]) %*% eta1[ij, , ] %*% tau_hat[j,])[1, 1]
    })}))
    
    #calcul de mu
    phi1[1] <- Psi1_hat %*% vec_S / sum(Psi1_hat)
    phi0[1] <- Psi0_hat %*% vec_S / sum(Psi0_hat)
    phi1[2] <- Psi1_hat %*% (vec_S^2) / sum(Psi1_hat) - phi1[1]^2
    phi0[2] <- Psi0_hat %*% (vec_S^2) / sum(Psi0_hat) - phi0[1]^2

    #calcul de fu:densité des 
    lfu <-cbind(dnorm(vec_S,phi0[1],sqrt(phi0[2]),log = T),dnorm(vec_S, phi1[1],sqrt(phi1[2]),log = T))
    dlfu = lfu[, 1] - lfu[, 2]
    # Troncature pour eviter les pb numeriques
    dlfu[which(abs(dlfu) > 100)] = sign(dlfu[which(abs(dlfu) > 100)]) * 100
    eta1 = 1 / (1 + exp(dlfu) %o% ((1 - Gamma_hat) / Gamma_hat))
    eta0 = 1 - eta1
    # Lissage des proba a conditionnelles
    eta0 = eta0 + epsilon_eta; eta1 = eta1 + epsilon_eta; 
    eta = eta0 + eta1; eta0 = eta0 / eta; eta1 = eta1 / eta
    
    #calcul de A
    A <-  eta0 * (rep(1, N)%o%log(1 - Gamma_hat)) + eta0 * lfu[,1] + 
      eta1 * (rep(1, N)%o%log(Gamma_hat)) + eta1 * lfu[,2]
    
    #calcul de la borne_inf pour control:
    vec_BI[(3*t)-1] =  borne_inf(tau_hat,Pi_hat,A,eta0,eta1)
    
    #########################################################################################################################
    #########################################################################################################################
    #etape VE
    #calcul de tau
    tau_old = tau_hat
    tau_new = matrix(0, n, K)
    for (j in 1:100){
      
      sapply(1:n, function(i){
        
        # m1.i = ensemble des paires dont i est le premier element
        m1.i = which(transfo_indices[, 1]==i); if(length(m1.i)>0){B1.i = array(A[m1.i, , ], dim=c(length(m1.i), K, K))}
        # m2.i = ensemble des paires dont i est le second element
        m2.i = which(transfo_indices[, 2]==i); 
        # Transpoition de A quand i est le second element de la paire m
        if(length(m2.i)>0){B2.i = array(A[m2.i, , ], dim=c(length(m2.i), K, K)); sapply(1:length(m2.i), function(m){B2.i[m, , ] <<- t(B2.i[m, , ])})}
        B.i = array(dim=c(length(m1.i)+length(m2.i), K, K))
        if(length(m1.i)>0){B.i[1:length(m1.i), , ] = B1.i}
        if(length(m2.i)>0){B.i[(length(m1.i)+1):(length(m1.i)+length(m2.i)), , ] = B2.i}
        #  Calcul de tau
        tau.j = tau_old[-i, ]; tau.i = rep(0, K)
        sapply(1:K, function(k){
          tau.i[k] <<- log(Pi_hat[k]) + sum(tau.j * B.i[, k, ])
        })
        tau.i = tau.i - max(tau.i); tau.i[which(tau.i < -100)] = -100
        tau.i = exp(tau.i); tau.i = tau.i / sum(tau.i)
        tau.i = tau.i + 1e-4; tau.i = tau.i / sum(tau.i)
        tau_new[i, ] <<- tau.i
      })
      tau_old = tau_new
    }
    
    # Test et mise a jour
    diff = max(abs(tau_new - tau_hat))
    if (t == (niter-1)){print("Maximum number of iterations reached")}
    tau_hat = tau_new
    
    #calcul de la borne_inf pour control:

    vec_BI[3*t] =  borne_inf(tau_hat,Pi_hat,A,eta0,eta1)
  # 
  }

  ############" reorder
  ord <- order(diag(Gamma_hat), decreasing  = TRUE)
  Pi_hat <- Pi_hat[ord]  
  output <- list(tau_init  = tau_init[,ord], tau  = tau_hat[,ord], 
                 phi0 = phi0, phi1 = phi1,
                 borne_inf = vec_BI[1:(3*t)])
  output$Pi_hat <- Pi_hat
  output$Gamma_hat <- Gamma_hat[ord,ord]
  output$Psi0 <- Psi0_hat
  output$Psi1 <- Psi1_hat
  
  return(output)
}


##########################################################################
# Calcul de la borne inf
borne_inf <- function(tau_hat,Pi_hat,A,eta0,eta1){
   borne_inf <- sum(tau_hat %*% log(Pi_hat))  - sum(tau_hat * log(tau_hat)) +
      .5*(sum(sapply(1:K, function(k){ sapply(1:K, function(l){
         t(tau_hat[,k]) %*% as.matrix(vect_mat_low(A[,k,l])) %*% tau_hat[,l]})
      }))) - 
      .5*(sum(sapply(1:K, function(k){ sapply(1:K, function(l){
         t(tau_hat[,k]) %*% 
            as.matrix(vect_mat_low(eta0[, k, l]*log(eta0[, k, l]) +
                                      eta1[, k, l]*log(eta1[, k, l]))) %*% 
            tau_hat[,l]})
      }))) 
   return(borne_inf)
}

