# PostDiff and Top-2 TS simulation

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

# df_ts is a matrix of PostDiff TS simulations. it's 41 by 6 by 5 by 3
# 41 corresponds to the tunning prarmeter 'c' goes from 0 to 1, with a increment of 0.025
# 6 corresponds to sample size being 100, 197, 300, 600, 785, and 1000, respectively
# 5 corresponds to effect size goes being 0.1, 0.2, 0.3, 0.5, and 0.8 (though 0.5 and 0.8 is not used in our final results)
# 3 corresponds to FPR (when calculating FPR, no matter which effect size dimension we are in, the actual effect size is set to 0), Power and reward.

# df_eps is a matrix of Top2 TS simulations. it's 41 by 6 by 5 by 3
# 41 corresponds to the tunning prarmeter 'beta' goes from 1 to 0.5, with a increment of -0.0125
# 6 corresponds to sample size being 100, 197, 300, 600, 785, and 1000, respectively
# 5 corresponds to effect size goes being 0.1, 0.2, 0.3, 0.5, and 0.8 (though 0.5 and 0.8 is not used in our final results)
# 3 corresponds to FPR (when calculating FPR, no matter which effect size dimension we are in, the actual effect size is set to 0), Power and reward.


