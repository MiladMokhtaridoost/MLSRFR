MLSRFR <- function(X1, X2, Y1, Y2, lam1, lam2, lam3, lam4, ranks)  {
  
  R1 <- ranks$R1_u
  R2 <- ranks$R2_u
  Rs <- ranks$Rs_u
  
  wthresh <- function(x, t){
    Wx <- abs(x) - t
    Wx <- (abs(Wx)+Wx)/2
    x <- sign(x) * Wx
    return(x)
  }
  maxIter <- 1000
  obj_vals <- matrix(0,nrow = maxIter, ncol=1)
  tol <- 0.0001
  D <- ncol(X1)
  T <- ncol(Y1)

  Wsx <- matrix(runif(Rs * D), ncol = Rs, nrow = D)
  Wsy <- matrix(runif(Rs * T), ncol = T, nrow = Rs)
  Wx1 <- matrix(runif((R1) * D), ncol =(R1), nrow = D) 
  Wx2 <- matrix(runif((R2) * D), ncol =(R2), nrow = D)
  Wy1 <- matrix(runif((R1) * T), ncol = T, nrow = (R1))
  Wy2 <- matrix(runif((R2) * T), ncol = T, nrow = (R2))
  W1sx <- cbind(Wsx, Wx1)
  W1sy <- rbind(Wsy, Wy1)
  W2sx <- cbind(Wsx, Wx2)
  W2sy <- rbind(Wsy, Wy2)
  
  obj_vals[1] <- 99999999999
  
  iter=1
  while (iter <= maxIter) {
    iter = iter+1 
    ########## updating Wy ######
    print(iter)

    Beta_1 <- norm((crossprod(Wx1,t(X1))) %*% (X1 %*% cbind((Wsx%*% Wsy), Wx1)),'F') + 2 * lam4
    Beta_2 <- norm((crossprod(Wx2,t(X2))) %*% (X2 %*% cbind((Wsx%*% Wsy), Wx2)),'F') + 2 * lam4
    Beta_s <- norm((crossprod(Wsx,t(X1))) %*% (X1 %*% cbind(Wsx, (Wx1%*% Wy1))),'F') + 
      norm((crossprod(Wsx,t(X2))) %*% (X2 %*% cbind(Wsx, (Wx2%*% Wy2))),'F') + 2 * lam4
    
    Ghat_1 <- (crossprod(Wx1,t(X1))) %*% ((X1 %*% W1sx) %*% W1sy - Y1) + 2 * lam4 * (Wy1)
    Ghat_2 <- (crossprod(Wx2,t(X2))) %*% ((X2 %*% W2sx) %*% W2sy - Y2) + 2 * lam4 * (Wy2)
    Ghat_s <- ((crossprod(Wsx,t(X1))) %*% ((X1 %*% W1sx) %*% W1sy - Y1)) +
      ((crossprod(Wsx,t(X2))) %*% ((X2 %*% W2sx) %*% W2sy - Y2)) + 2 * lam4 * (Wsy)
    Wy1 <- wthresh(Wy1 - Ghat_1 / Beta_1, lam3 / Beta_1)
    Wy2 <- wthresh(Wy2 - Ghat_2 / Beta_2, lam3 / Beta_2)
    Wsy <- wthresh(Wsy - Ghat_s / Beta_s, lam3 / Beta_s)

    W1sy <- rbind(Wsy, Wy1)
    W2sy <- rbind(Wsy, Wy2)
    
    ########### updating Wx ############
    alpha_1 <- norm(crossprod(X1) %*% Wsx, 'F') * norm(tcrossprod(W1sy, Wy1), 'F') + 2 * lam2
    alpha_2 <- norm(crossprod(X2) %*% Wsx, 'F') * norm(tcrossprod(W2sy, Wy2), 'F') + 2 * lam2
    alpha_s <- norm(crossprod(X1) %*% Wx1, 'F') * norm(tcrossprod(W1sy, Wsy), 'F') +
      norm(crossprod(X2) %*% Wx2, 'F') * norm(tcrossprod(W2sy, Wsy), 'F') + 2 * lam2
    
    Hhat_1 <- -t(X1) %*% (Y1 %*% t(Wy1)) + t(X1) %*% ((X1 %*% W1sx) %*% (tcrossprod(W1sy ,Wy1))) + 2 * lam2 * (Wx1) 
    Hhat_2 <- -t(X2) %*% (Y2 %*% t(Wy2)) + t(X2) %*% ((X2 %*% W2sx) %*% (tcrossprod(W2sy ,Wy2))) + 2 * lam2 * (Wx2)
    Hhat_s <- -t(X1) %*% (Y1 %*% t(Wsy)) + t(X1) %*% ((X1 %*% W1sx) %*% (tcrossprod(W1sy ,Wsy))) +
      -t(X2) %*% (Y2 %*% t(Wsy)) + t(X2) %*% ((X2 %*% W2sx) %*% (tcrossprod(W2sy ,Wsy))) + 2 * lam2 * (Wsx)
    
    Wx1 <- wthresh(Wx1 - Hhat_1 / alpha_1, lam1 / alpha_1)
    Wx2 <- wthresh(Wx2 - Hhat_2 / alpha_2, lam1 / alpha_2)
    Wsx <- wthresh(Wsx - Hhat_s / alpha_s, lam1 / alpha_s) 

    W1sx <- cbind(Wsx, Wx1)
    W2sx <- cbind(Wsx, Wx2)
    
    ########### stopping check ###########
    obj_vals[iter] <- 0.5 * norm(Y1 - (X1 %*% W1sx) %*% W1sy, 'F')^2 + norm(Y2 - (X2 %*% W2sx) %*% W2sy, 'F')^2 + 
      lam1 * (sum(abs(Wx1)) + sum(abs(Wx2)) + sum(abs(Wsx))) + lam2 * (sqrt(sum(Wx1^2)) + sqrt(sum(Wx2^2)) + sqrt(sum(Wsx^2))) + 
      lam3 * (sum(abs(Wy1)) + sum(abs(Wy2)) + sum(abs(Wsy))) + lam4 * (sqrt(sum(Wy1^2)) + sqrt(sum(Wy2^2)) + sqrt(sum(Wsy^2))) 
     
    
    dif <- (abs(obj_vals[iter]- obj_vals[iter-1]))/obj_vals[iter-1] 
    if (dif < tol || is.nan(dif)) {
      if (((qr(Wy1)$rank != R1) || (qr(Wx1)$rank != R1))) {
        R1 <- R1-1
        
        Wx1 <- matrix(runif((R1) * D), ncol =(R1), nrow = D) 
        Wy1 <- matrix(runif((R1) * T), ncol = T, nrow = (R1))
        W1sx <- cbind(Wsx, Wx1)
        W1sy <- rbind(Wsy, Wy1)
        
        if (((qr(Wy2)$rank != R2) || (qr(Wx2)$rank != R2))) {
          R2 <- R2-1
          
          Wx2 <- matrix(runif((R2) * D), ncol =(R2), nrow = D)
          Wy2 <- matrix(runif((R2) * T), ncol = T, nrow = (R2))
          W2sx <- cbind(Wsx, Wx2)
          W2sy <- rbind(Wsy, Wy2)
          
          if (((qr(Wsy)$rank != Rs) || (qr(Wsx)$rank != Rs))) {
            Rs <- Rs-1
            
            Wsx <- matrix(runif(Rs * D), ncol = Rs, nrow = D)
            Wsy <- matrix(runif(Rs * T), ncol = T, nrow = Rs)
            W1sx <- cbind(Wsx, Wx1)
            W1sy <- rbind(Wsy, Wy1)
            W2sx <- cbind(Wsx, Wx2)
            W2sy <- rbind(Wsy, Wy2)
          }
        }
      }
      else  {
        
        break
      }
    }
  }
  rank_star <- list()
  rank_star$R1 <- R1
  rank_star$R2 <- R2
  rank_star$Rs <- Rs

  return(list(Wx1_best = Wx1, Wx2_best = Wx2, Wsx_best = Wsx, Wy1_best = Wy1, Wy2_best = Wy2, Wsy_best = Wsy, R_stars = rank_star, iteration=iter, obj_vals=obj_vals))
}
