

######################################################### Define PostDiff TS  and Top-2 TS
TSPostDiff_sim <- function(pa,pb,n,c){
  # n is the sample size (785 or 197 in our case)
  arm_a_successes <- 0
  arm_a_failures <- 0
  arm_b_successes <- 0
  arm_b_failures <- 0
  na <- c()
  sa <- c()
  sb <- c()
  for (i in 1:n){
    draw_a <- rbeta(1,1+arm_a_successes,1+arm_a_failures)
    draw_b <- rbeta(1,1+arm_b_successes,1+arm_b_failures)
    if (abs(draw_a-draw_b)<c){
      draw_a <- runif(1)
      draw_b <- runif(1)
    } else {
      draw_a <- rbeta(1,1+arm_a_successes,1+arm_a_failures)
      draw_b <- rbeta(1,1+arm_b_successes,1+arm_b_failures)
    }
    if (draw_a>draw_b){
      if(runif(1)<pa){
        arm_a_successes <- arm_a_successes+1
      } else {
        arm_a_failures <- arm_a_failures+1
      }
    } else{
      if(runif(1)<pb){
        arm_b_successes <- arm_b_successes+1
      } else {
        arm_b_failures <- arm_b_failures+1
      }
    }
    na[i] <- arm_a_successes+arm_a_failures
    sa[i] <- arm_a_successes
    sb[i] <- arm_b_successes
  }
  #na: number of times arm A was pulled. nb: same
  
  
  return(list(na=na,nb=c(1:n)-na,sa=sa,sb=sb))
  #return(list(WaldScore=WaldScore,na=na,nb=nb,sa=arm_a_successes,sb=arm_b_successes))
}


Top2_TS_sim <- function(pa,pb,n,epsilon){
  # n is the sample size (785 or 197 in our case)
  arm_a_successes <- 0
  arm_a_failures <- 0
  arm_b_successes <- 0
  arm_b_failures <- 0
  na <- c()
  sa <- c()
  sb <- c()
  for (i in 1:n){
    draw_a <- rbeta(1,1+arm_a_successes,1+arm_a_failures)
    draw_b <- rbeta(1,1+arm_b_successes,1+arm_b_failures)
    if (runif(1)<epsilon){
      draw_a <- runif(1)
      draw_b <- runif(1)
    }
    
    if (draw_a>draw_b){
      if(runif(1)<pa){
        arm_a_successes <- arm_a_successes+1
      } else {
        arm_a_failures <- arm_a_failures+1
      }
    } else{
      if(runif(1)<pb){
        arm_b_successes <- arm_b_successes+1
      } else {
        arm_b_failures <- arm_b_failures+1
      }
    }
    na[i] <- arm_a_successes+arm_a_failures
    sa[i] <- arm_a_successes
    sb[i] <- arm_b_successes
  }
  #na: number of times arm A was pulled. nb: same
  return(list(na=na,nb=c(1:n)-na,sa=sa,sb=sb))
}




##################################################### Simulation start
lh <- c(100,197,300,600,785,1000)
h <- length(lh)
n <- max(lh)
m<-41
e_s <- c(0.1,0.2,0.3,0.5,0.8)
hh <- length(e_s)
start_time = Sys.time()
B <- 10000
df_eps <- array(0,dim=c(m,h,hh,3))
df_ts <- array(0,dim=c(m,h,hh,3))



kk <- 0
for (s in 1:hh){
  k <- kk
  for (j in 1:B){
    effect_size=e_s[s]
    
    set.seed(k+10000)
    df_eps0 <- array(dim=c(m,h,3))
    df_ts0 <- array(dim=c(m,h,3))
    for (i in 1:m){
      x <- (i-1)/(m-1)
      # FPR
      r1 <- Top2_TS_sim(0.5,0.5,n,x)
      r2 <- TSPostDiff_sim(0.5,0.5,n,x)
      est_a1 <- r1$sa/r1$na
      est_b1 <- r1$sb/r1$nb
      
      est_a2 <- r2$sa/r2$na
      est_b2 <- r2$sb/r2$nb
      w_eps <- (est_a1-est_b1)/sqrt(est_a1*(1-est_a1)/(r1$na) + est_b1*(1-est_b1)/(r1$nb))
      w_ts <- (est_a2-est_b2)/sqrt(est_a2*(1-est_a2)/(r2$na) + est_b2*(1-est_b2)/(r2$nb))
      
      df_eps0[i,,1] <- (abs(w_eps)<1.96)[lh]
      df_ts0[i,,1] <- (abs(w_ts)<1.96)[lh]
      
      #Power, Reward
      r1 <- Top2_TS_sim(0.5+effect_size/2,0.5-effect_size/2,n,x)
      r2 <- TSPostDiff_sim(0.5+effect_size/2,0.5-effect_size/2,n,x)
      est_a1 <- r1$sa/r1$na
      est_b1 <- r1$sb/r1$nb
      
      est_a2 <- r2$sa/r2$na
      est_b2 <- r2$sb/r2$nb
      w_eps <- (est_a1-est_b1)/sqrt(est_a1*(1-est_a1)/(r1$na) + est_b1*(1-est_b1)/(r1$nb))
      w_ts <- (est_a2-est_b2)/sqrt(est_a2*(1-est_a2)/(r2$na) + est_b2*(1-est_b2)/(r2$nb))
      
      r_eps <- 0.5-effect_size/2+r1$na*effect_size/c(1:n)
      r_ts <- 0.5-effect_size/2+r2$na*effect_size/c(1:n)
      
      df_eps0[i,,2] <- (w_eps>1.96)[lh]
      df_ts0[i,,2] <- (w_ts>1.96)[lh]
      
      df_eps0[i,,3] <- r_eps[lh]
      df_ts0[i,,3] <- r_ts[lh]
      df_eps0[is.na(df_eps0)] <- 0
      df_ts0[is.na(df_ts0)] <- 0
    }  
    df_eps[,,s,] <- (df_eps[,,s,]*k+df_eps0)/(k+1)
    df_ts[,,s,] <- (df_ts[,,s,]*k+df_ts0)/(k+1)
    k <- k+1
  }
}
##################################################### Simulation end




############################################################################# 9-2 FPR unifomity
# update 9-5
df_eps3 <- df_eps0
df_ts3 <- df_ts0
set.seed(1)
m <- 41

B <- 10000
df_eps0 <- array(dim=c(41,B,6,2))
df_ts0 <- array(dim=c(41,B,6,2))

df_eps0[c(1:21)*2-1,,,] <- df_eps3[,1:B,,]
df_ts0[c(1:21)*2-1,,,] <- df_ts3[,1:B,,]


set.seed(1)
n <- 1000
lh <- c(100,197,300,600,785,1000)
for (j in 1:10000){
  for (i in 1:m){
    x <- (i-1)/(m-1)
    # FPR
    r1 <- Top2_TS_sim(0.5,0.5,n,x)
    r2 <- TSPostDiff_sim(0.5,0.5,n,x)
    
    est_a1 <- r1$sa/r1$na
    est_b1 <- r1$sb/r1$nb
    
    est_a2 <- r2$sa/r2$na
    est_b2 <- r2$sb/r2$nb
    
    w_eps <- (est_a1-est_b1)/sqrt(est_a1*(1-est_a1)/(r1$na) + est_b1*(1-est_b1)/(r1$nb))
    w_ts <- (est_a2-est_b2)/sqrt(est_a2*(1-est_a2)/(r2$na) + est_b2*(1-est_b2)/(r2$nb))
    
    df_eps0[i,j,,1] <- (abs(w_eps)<1.96)[lh]
    df_ts0[i,j,,1] <- (abs(w_ts)<1.96)[lh]
    
    df_eps0[i,j,,2] <- r1$na[lh]/lh
    df_ts0[i,j,,2] <- r2$na[lh]/lh
    df_eps0[i,j,,3] <- est_a1>est_b1
    df_ts0[i,j,,3] <- est_a2>est_b2
    
    
  }  
  
}


eps0 <- df_eps0[,,1,2]
ts0 <- df_ts0[,,1,2]
ts0[ts0>0.5] <- 1-ts0[ts0>0.5]
eps0[eps0>0.5] <- 1-eps0[eps0>0.5]

apply(ts0,1,mean)[tsl]
apply(eps0,1,mean)[epl]

plot(density(ts0[3,]))

cvalue <- 3
evalue <- 5
samplesize <- 4
p1 <- hist(df_ts0[cvalue,,samplesize,2])                     # centered at 4
p2 <- hist(df_eps0[evalue,,samplesize,2])                     # centered at 6
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1),main=paste0('Sample Size = ', c(100,197,300,600,785,1000)[samplesize]),xlab = 'Propotion that Arm A is Sampled')  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T) 
legend("topright", c(paste0("TS PostDiff, c=0.",(cvalue-1)*5 ), paste0("epsilon TS, eps=0.",(evalue-1)*5)), col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), lwd=10)


eps0 <- df_eps0[,,2]
ts0 <- df_ts0[,,2]
sqrt(apply(ts0,1,var)[c(3,5,7,9)]/10000)

sqrt(apply(ts0,1,var)[c(3,5,7,9)]/10000)
apply(eps0,1,var)[c(5,9,13,17)]

sqrt(apply(ts0,1,var)[c(3,5,7,9)]*(1-apply(ts0,1,var)[c(3,5,7,9)])/50000)


#var of FPR
apply(ts0,2,mean)[tsl]

eps1 <- df_eps0[,,4,1]
ts1 <- df_ts0[,,4,1]



1-apply(eps1,1,mean)[epl]

sqrt(apply(eps1,1,mean)[epl]*(1-apply(eps1,1,mean)[epl])/50000)

1-apply(ts1,1,mean)[tsl]
sqrt((1-apply(ts1,1,mean)[tsl])*apply(ts1,1,mean)[tsl]/50000)

for (j in 1:B){
  for (i in 1:m){
    x <- (i-1)/(m-1)
    # FPR
    r1 <- Top2_TS_sim(0.5,0.5,n,x)
    r2 <- TSPostDiff_sim(0.5,0.5,n,x)
    
    est_a1 <- r1$sa/r1$na
    est_b1 <- r1$sb/r1$nb
    
    est_a2 <- r2$sa/r2$na
    est_b2 <- r2$sb/r2$nb
    
    w_eps <- (est_a1-est_b1)/sqrt(est_a1*(1-est_a1)/(r1$na) + est_b1*(1-est_b1)/(r1$nb))
    w_ts <- (est_a2-est_b2)/sqrt(est_a2*(1-est_a2)/(r2$na) + est_b2*(1-est_b2)/(r2$nb))
    
    df_eps0[i,j,1] <- (abs(w_eps)<1.96)[n]
    df_ts0[i,j,1] <- (abs(w_ts)<1.96)[n]
    
    df_eps0[i,j,2] <- r1$na[n]/n
    df_ts0[i,j,2] <- r2$na[n]/n
    
    
  }  
  
}

a1 <- apply(df_eps0[,,1],1,mean)

a2 <- apply(df_ts0[,,1],1,mean)
a3 <- apply(df_eps0[,,2],1,var)
a4 <- apply(df_ts0[,,2],1,var)
b1 <- cbind(a1,a3,a2,a4)
write.csv(b1,file='b1.csv')

plot(apply(df_eps0,c(1,3),mean))
plot(apply(df_ts0,c(1,3),mean))

plot(density(df_eps0[1,,2]))
for ( i in 1:11){
  plot(density(df_ts0[i,,2]))
}


###############################################################

end_time = Sys.time()
end_time-start_time
save(df_eps,df_ts,k,file = 'TSsim1.RData')
kk <- k

apply(sqrt((df_ts[,,1])^2+(df_ts[,,2])^2+(df_ts[,,3])^2),2,sum)

apply(sqrt((df_eps[,,1])^2+(df_eps[,,2])^2+(df_eps[,,3])^2),2,sum)

write.csv((df_ts[,,]),file='ts.csv')
write.csv((df_eps[,,]),file='eps.csv')



ef <- 1
sam <- 1

lh <- c(100,197,300,600,785,1000)
h <- length(lh)
n <- max(lh)
m<-41
e_s <- c(0.1,0.2,0.3,0.5,0.8)

df_eps[,sam,ef,]
df_ts[,sam,ef,]

df_eps[,sam,,1]
apply(df_eps[,sam,,1],1,mean)
em <- apply(df_eps[,,,1],c(1,2),mean)
emt <- apply(df_ts[,,,1],c(1,2),mean)
df_eps1 <- df_eps
df_ts1 <- df_ts

for (i in 1:5){
  df_eps[,,i,1] <- em
  df_ts[,,i,1] <- emt
}

df_eps[,1,,1]
df_ts[,1,,1]


df_eps[,,,1] <- 1-df_eps[,,,1]
df_ts[,,,1] <- 1-df_ts[,,,1]

##########################
fpr1 <- 6
fpr <- fpr1/100
sam <- 1
ef <- 1

for ( sam in 1:6){
  for (ef in 1:5){
    for (fpr1 in c(6)){
      k <- df_eps[df_eps[,sam,ef,1]<fpr,sam,ef,c(2,3)]
      kk <- df_ts[df_ts[,sam,ef,1]<fpr,sam,ef,c(2,3)]
      xmin <- min(min(k[,1],min(kk[,1])))
      xmax <- max(max(k[,1],max(kk[,1])))
      ymin <- min(min(k[,2],min(kk[,2])))
      ymax <- max(max(k[,2],max(kk[,2])))
      ep <- 0.005
      
      mt <- paste0(' Sample Size = ' ,c(100,197,300,600,785,1000)[sam], ', Effect Size = ', c(0.1,0.2,0.3,0.5,0.8)[ef],', FPR bar = ', fpr1, '%' )
      df <- array(k,dim=dim(k))
      dff <- array(kk,dim=dim(kk))
      jpeg(paste0(sam,ef,fpr1,'.jpg'))
      
      plot(df[,1], df[,2], main = mt,
           xlab = "Power", ylab = "Reward",
           xlim = c(xmin-ep/2,xmax+ep), ylim=c(ymin-ep/2,ymax+ep/2),
           pch = 19,  col='blue')
      
      points(dff[,1], dff[,2],pch = 19 , col='red')
      legend(xmax-ep/2,ymax, legend=c("TSPostDiff", "Epsilon Greedy"),
             col=c("red", "blue"),pch=c(19,19), cex=0.8)
      dev.off()
    }
  }
}









for ( sam in 1:6){
  for (ef in 1:5){
    for (fpr1 in c(8)){
      fpr <- fpr1/100
      l <- c(1:11)*4-3
      k <- cbind(df_eps[l,sam,ef,],c(0:10)/10)
      kk <- cbind(df_ts[l,sam,ef,],c(0:10)/10)
      
      k <- k[k[,1]<fpr,]
      kk <- kk[kk[,1]<fpr,]
      
      
      df <- array(k,dim=dim(k))
      df <- data.frame(df)
      colnames(df) <- c('fpr','power','reward','c')
      dff <- array(kk,dim=dim(kk))
      dff <- data.frame(dff)
      colnames(dff) <- c('fpr','power','reward','c')
      xmin <- min(min(k[,2],min(kk[,2])))
      xmax <- max(max(k[,2],max(kk[,2])))
      ymin <- min(min(k[,3],min(kk[,3])))
      ymax <- max(max(k[,3],max(kk[,3])))
      mt <- paste0(' Sample Size = ' ,c(100,197,300,600,785,1000)[sam], ', Effect Size = ', c(0.1,0.2,0.3,0.5,0.8)[ef],', FPR bar = ', fpr1, '%' )
      
      
      df_kep <- c()
      for ( i in 1:dim(df)[1]){
        s <- sum((df$power>=df$power[i]) & (df$reward>=df$reward[i]))+sum((dff$power>=df$power[i]) & (dff$reward>=df$reward[i]))
        df_kep[i] <- (s==1)
      }
      
      dff_kep <- c()
      for ( i in 1:dim(dff)[1]){
        s <- sum((df$power>=dff$power[i]) & (df$reward>=dff$reward[i]))+sum((dff$power>=dff$power[i]) & (dff$reward>=dff$reward[i]))
        dff_kep[i] <- (s==1)
      }
      
      df <- df[df_kep,]
      # dff is ts
      dff <- dff[dff_kep,]
      jpeg(paste0(sam,ef,fpr1,'.jpg'))
      plot(dff[,2], dff[,3], main = mt,
           xlab = "Power", ylab = "Reward",
           xlim = c(xmin-ep/2,xmax+ep), ylim=c(ymin-ep/2,ymax+ep/2),
           pch = 19,  col='red1',cex=3)
      if(dim(df)[1]>0){
        points(df[,2], df[,3],pch = 19 , col='lightblue',cex=3)
      }
      
      #legend(xmax-ep/2,ymax, legend=c("TSPostDiff", "Epsilon Greedy"),
      #       col=c("red", "blue"),pch=c(19,19), cex=0.8)
      
      if(dim(dff)[1]>0){
        text(reward~power, labels=c,data=dff, cex=1.5, font=1)
      }
      if(dim(df)[1]>0){
        text(reward~power, labels=c,data=df, cex=1.5, font=1)
      }
      
      dev.off()
      
    }
  }
}















































for ( sam in 1:6){
  for (ef in 1:5){
    for (fpr1 in c(8)){
      fpr <- fpr1/100
      #l <- c(1:11)*4-3
      #k <- cbind(df_eps[l,sam,ef,],c(0:10)/10)
      #kk <- cbind(df_ts[l,sam,ef,],c(0:10)/10)
      
      
      
      k <- cbind(c(0:40)/40, df_eps[,sam,ef,],rep('',41),c(0:40)/40, df_ts[,sam,ef,])
      
      #k <- k[k[,1]<fpr,]
      #kk <- kk[kk[,1]<fpr,]
      
      
      df <- array(k,dim=dim(k))
      df <- data.frame(df)
      colnames(df) <- c('epsilon','FPR','Power','Reward','','c','FPR','Power','Reward')
      mt <- paste0(' Sample Size = ' ,c(100,197,300,600,785,1000)[sam], ', Effect Size = ', c(1,2,3,5,8)[ef] )
      write.csv(df,file=paste0(mt,'.csv'))
      
      
    }
  }
}


library(plotly)



plot_ly(x=data$x, y=data$y, z=data$z, type="scatter3d", mode="markers", color=c(1:21)/21)


for ( sam in 1:6){
  for (ef in 1:5){
    for (fpr1 in c(8)){
      fpr <- fpr1/100
      #l <- c(1:11)*4-3
      #k <- cbind(df_eps[l,sam,ef,],c(0:10)/10)
      #kk <- cbind(df_ts[l,sam,ef,],c(0:10)/10)
      
      
      
      k <- cbind(c(0:40)/40, df_eps[,sam,ef,],rep(1,41))
      kk <- cbind(c(0:40)/40, df_ts[,sam,ef,],rep(2,41))
      #k <- k[k[,1]<fpr,]
      #kk <- kk[kk[,1]<fpr,]
      k1 <- rbind(k,kk)
      
      df <- array(k1,dim=dim(k1))
      df <- data.frame(df)
      df$X2 <- 1-df$X2
      colnames(df) <- c('c','1 - FPR','Power','Reward','color')
      
      jpeg(paste0('3d',sam,ef,'.jpg'))
      plot_ly(x=df$`1 - FPR`, y=df$Power, z=df$Reward, type="scatter3d", mode="markers", color=df$color)
      dev.off()
    }
  }
}




