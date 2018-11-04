spls.coeff<- function(tf.exp,pw.exp) {
  nTF <- dim(tf.exp)[2]
  nPW <- dim(pw.exp)[2]
  sample_size<-dim(tf.exp)[1]

  cv <-cv.spls(tf.exp,pw.exp,eta = seq(0.6,0.9,0.1),K = c(5:10), 
               kappa=0.5, select="pls2", fit="simpls",scale.x=TRUE, scale.y=TRUE, plot.it=F)
  eta = cv$eta.opt
  K = cv$K.opt
  f.ori <- spls(tf.exp,pw.exp,eta = eta,K = K)
  ci.f <- ci.spls(f.ori,coverage = 0.95,plot.it = F,plot.var = F)
  f.cor.mv <- correct.spls(ci.f,plot.it = F)
  coef.f.corrected.mv <- abs(f.cor.mv)
  sorted_score.cor.mv = data.frame(matrix(nrow = nTF,ncol = 0))#make this dynamic
  pw_genes <-colnames(pw.exp)
 
  for (pw_i in 1:length(pw_genes)){
    temp<-data.frame(rownames(coef.f.corrected.mv),coef.f.corrected.mv[,pw_i])
    colnames(temp)<-c(paste(pw_genes[pw_i],"pwg",sep = "_"),paste(pw_genes[pw_i],"coeff",sep = "_"))
    temp<-temp[order(-temp[2]),]
    sorted_score.cor.mv <- cbind(sorted_score.cor.mv,temp)
    
  }
  
  return(sorted_score.cor.mv)
  
}


coeff.boot<-function(tf.exp,pw.exp){
  nTF <- dim(tf.exp)[2]
  nPW <- dim(pw.exp)[2]
  sample_size<-dim(tf.exp)[1]
  cv <-cv.spls(tf.exp,pw.exp,eta = seq(0.6,0.9,0.1),K = c(5:10), 
               kappa=0.5, select="pls2", fit="simpls",scale.x=TRUE, scale.y=TRUE, plot.it=F)
  eta = cv$eta.opt
  K = cv$K.opt
  iters = 1e3
  f.ori <- spls(tf.exp,pw.exp,eta = eta,K = K)
  ci.f <- ci.spls(f.ori,coverage = 0.95,plot.it = F,plot.var = F)
  f.cor.mv <- correct.spls(ci.f,plot.it = F)
  coef.f.corrected.mv <- abs(f.cor.mv)
  
  cl<-makeCluster(8)
  registerDoParallel(cl)
  
  ls<-foreach(icount(iters),.packages = "spls") %dopar% {
    boot_ind<-sample(1:sample_size,replace = T)
    boot_tf.exp<-tf.exp[boot_ind,]
    boot_pw.exp<-pw.exp[boot_ind,]
    f.boot <- spls(boot_tf.exp,boot_pw.exp,eta = eta,K = K)
    ci.f <- ci.spls(f.boot,coverage = 0.95,plot.it = F,plot.var = F)
    cf <- correct.spls(ci.f,plot.it = F)
    coef.f.corrected.mv <- abs(cf)
    b<-as.matrix(coef.f.corrected.mv)
    b
  }
  cc <- matrix(0, ncol=nPW, nrow=nTF, byrow=T)
  
  
  for(i in 1:dim(cc)[1]){
    for(j in 1:dim(cc)[2]){
      betas<-c()
      for(k in 1:iters){
        betas<-c(betas,ls[[k]][i,j])
      }
      beta <- coef.f.corrected.mv[i,j]
      #bootSigma <- sqrt(sum((betas-bootmean)^2)/iters)
      booter <- sd(betas)/sqrt(128)
      
      if(booter!=0){
        
        tval <- beta/booter
        #bVar <-bootvar
        #pval<-2*pt(-abs(tval),df=127)
      }else{
        tval <-0
        #pval<-2*pt(-abs(tval),df=127)
      }
      cc[i,j] <- tval
      
    }
  }
  
  colnames(cc)<-colnames(pw.exp)
  rownames(cc)<-colnames(tf.exp)
  stopCluster(cl)
  
  sorted_score.cor.mv = data.frame(matrix(nrow = nTF,ncol = 0))#make this dynamic
  pw_genes <-colnames(pw.exp)
  
  for (pw_i in 1:length(pw_genes)){
    temp<-data.frame(rownames(cc),cc[,pw_i])
    colnames(temp)<-c(paste(pw_genes[pw_i],"pwg",sep = "_"),paste(pw_genes[pw_i],"t.statistic",sep = "_"))
    temp<-temp[order(-temp[2]),]
    sorted_score.cor.mv <- cbind(sorted_score.cor.mv,temp)
    
  }
  
  return(sorted_score.cor.mv)
  
}

coeff.perm.int <-function(tf.exp,pw.exp){
  nTF <- dim(tf.exp)[2]
  nPW <- dim(pw.exp)[2]
  sample_size<-dim(tf.exp)[1]
  
  cv <-cv.spls(tf.exp,pw.exp,eta = seq(0.8,0.9,0.1),K = c(8:10), 
               kappa=0.5, select="pls2", fit="simpls",scale.x=TRUE, scale.y=TRUE, plot.it=F)
  eta = cv$eta.opt
  K = cv$K.opt
  
  cc <- matrix(0, ncol = nPW, nrow=nTF, byrow=T)
  iters<-1e3
  cl<-makeCluster(8)
  registerDoParallel(cl)
  
  ls<-foreach(icount(iters),.packages = "spls") %dopar% {
    pw.exp.noise <- pw.exp[sample(nrow(pw.exp)),sample(ncol(pw.exp))]
    f.noise <- spls(tf.exp,pw.exp.noise,eta = eta,K = K)
    coef.f.noise <- abs(coef(f.noise))
    boolf.noise<-coef.f.noise > 0
    boolf.noise<-boolf.noise*1
    sparse.coeff.noise<-coef.f.noise*boolf.noise
    b<-as.matrix(sparse.coeff.noise)
    b
  }
  
  
  stopCluster(cl)
  
  for(i in 1:dim(cc)[1]){
    for(j in 1:dim(cc)[2]){
      betas<-c()
      for(k in 1:iters){
        betas<-c(betas,ls[[k]][i,j])
      }

      cc[i,j] <- mean(betas)
      
    }
  }
  
  colnames(cc)<-colnames(pw.exp)
  rownames(cc)<-colnames(tf.exp)
  
  sorted_score.cor.mv = data.frame(matrix(nrow = nTF,ncol = 0))#make this dynamic
  pw_genes <-colnames(pw.exp)
  
  for (pw_i in 1:length(pw_genes)){
    temp<-data.frame(rownames(cc),cc[,pw_i])
    colnames(temp)<-c(paste(pw_genes[pw_i],"pwg",sep = "_"),paste(pw_genes[pw_i],"t.statistic",sep = "_"))
    temp<-temp[order(-temp[2]),]
    sorted_score.cor.mv <- cbind(sorted_score.cor.mv,temp)
  }
  return(sorted_score.cor.mv)

}