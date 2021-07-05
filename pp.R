library(tidyverse)
# TODO: 2 figures
# 1: spaghetti plots by patient
# 2. diffs/ most likely graph

set.seed(1)
CD3s = read.csv("posterior.samples/CD3samples.csv")[,-1]
CD8s = read.csv("posterior.samples/CD8samples.csv")[,-1]
CD4s = read.csv("posterior.samples/CD4samples.csv")[,-1]
VLs = read.csv("posterior.samples/VLOADsamples.csv")[,-1]
Ys = read.csv("posterior.samples/Ysamples.csv")[,-1]

names(CD3s)[5:8] = names(CD8s)[5:8] = names(CD4s)[5:8] = names(VLs)[5:8] = c("b_0", "b_1", "b_2","sigsq") # b_2 is treatment

# y-variables: cbind(1, log(dat$CD4N_1), log(dat$CD8N_1), log(dat$CD3N_1), dat$VL_1, dat$A2_0, dat$A2_1),

dat = read.csv("combined2.csv")[,-1]
n = nrow(dat)

dat$VL_0 = ifelse(dat$VLOAD_0 < exp(10), 0, 1)
dat$VL_1 = ifelse(dat$VLOAD_1 < exp(10), 0, 1)

dat$y11 = NA
dat$y00 = NA
dat$t1 = ifelse(dat$A2_0 == 1 & dat$A2_1 == 1, 1, 0)
dat$t0 = ifelse(dat$A2_0 == 0 & dat$A2_1 == 0, 1, 0)

#
par(mfrow = c(2,1))
hist(filter(dat, t1 == 0)$CD4N_2, xlim = c(0,1500))
hist(filter(dat, t1 == 1)$CD4N_2, xlim = c(0, 1500))

# draws of pp
row = 201
ndraws = 50

diffs = array(0, dim = ndraws)
draws11 = array(0, dim = c(nrow(dat), ndraws))
draws00 = array(0, dim = c(nrow(dat), ndraws))

for(draw in 1:ndraws) {
  
  # CD3
  
  if ( CD3s[row+draw, 4] == 1 ) { # if there is an edge from treatment to CD3, impute CD3s[row+draw, 4] == 1
    
    cd31 = cbind(1, log(dat$CD3N_0), 1) %*% as.numeric( CD3s[row+draw, 5:7])  # XB (treated at time 1)
    cd3draw1 = rnorm(n, mean = cd31, sd = sqrt( CD3s[row + draw, 8]) ) # fixed treatment to 1
    cd3draw1 = ifelse(dat$A2_0 == 1, log(dat$CD3N_1), cd3draw1)
    
    
    cd30 = cbind(1, log(dat$CD3N_0), 0) %*% as.numeric( CD3s[row+draw, 5:7])  # XB (treated at time 1)
    print(head(cd30))
    cd3draw0 = rnorm(n, mean = cd30, sd = sqrt( CD3s[row + draw, 8]) ) # fixed treatment to 1
    cd3draw0 = ifelse(dat$A2_0 == 0, log(dat$CD3N_1), cd3draw0)
  } else {
    cd3draw1 = cd3draw0 = log(dat$CD3N_1)
  }
  
  # CD8
  if ( CD8s[row+draw, 4] == 1 ) { # if there is an edge from treatment to CD8, impute
    
    cd81 = cbind(1, log(dat$CD8N_0), 1) %*% as.numeric( CD3s[row+draw, 5:7])  # XB (treated at time 1)
    cd8draw1 = rnorm(n, mean = cd81, sd = sqrt( CD8s[row + draw, 8]) ) # fixed treatment to 1
    cd8draw1 = ifelse(dat$A2_0 == 1, log(dat$CD8N_1), cd8draw1)
    
    
    cd80 = cbind(1, log(dat$CD8N_0), 0) %*% as.numeric( CD8s[row+draw, 5:7])  # XB (treated at time 1)
    cd8draw0 = rnorm(n, mean = cd80, sd = sqrt( CD8s[row + draw, 8]) ) # fixed treatment to 1
    cd8draw0 = ifelse(dat$A2_0 == 0, log(dat$CD8N_1), cd8draw0)
  } else {
    cd8draw1 = cd8draw0 = log(dat$CD8N_1)
  }
  
  # CD4
  
  if ( CD4s[row+draw, 4] == 1 ) { # if there is an edge from treatment to CD8, impute
    
    cd41 = cbind(1, log(dat$CD4N_0), 1) %*% as.numeric( CD4s[row+draw, 5:7])  # XB (treated at time 1)
    cd4draw1 = rnorm(n, mean = cd41, sd = sqrt( CD4s[row + draw, 8]) ) # fixed treatment to 1
    cd4draw1 = ifelse(dat$A2_0 == 1, log(dat$CD4N_1), cd4draw1)
    
    
    cd40 = cbind(1, log(dat$CD4N_0), 0) %*% as.numeric( CD4s[row+draw, 5:7])  # XB (treated at time 1)
    cd4draw0 = rnorm(n, mean = cd40, sd = sqrt( CD4s[row + draw, 8]) ) # fixed treatment to 1
    cd4draw0 = ifelse(dat$A2_0 == 0, log(dat$CD4N_1), cd4draw0)
  } else {
    cd4draw1 = cd4draw0 = log(dat$CD4N_1)
  }
  
  
  # VL
  if ( VLs[row+draw, 4] == 1 ) { # if there is an edge from treatment to VL, impute
    
    VL1 = cbind(1, dat$VL_0, 1) %*% as.numeric(VLs[row+draw, 5:7])  # eta
    VLdraw1 = rbinom(n, 1, prob = plogis(VL1)) # fixed treatment to 1
    VLdraw1 = ifelse(dat$A2_0 == 1, dat$VL_0, VLdraw1)
    
    
    VL0 = cbind(1, dat$VL_1, 0) %*% as.numeric(VLs[row+draw, 5:7]) # eta
    VLdraw0 = rbinom(n, 1, prob = plogis(VL0)) # fixed treatment to 1
    VLdraw0 = ifelse(dat$A2_0 == 0, dat$VL_0, VLdraw0)
  } else {
    VLdraw1 = VLdraw0 = dat$VL_1
  }
  
  # y-variables: cbind(1, log(dat$CD4N_1), log(dat$CD8N_1), log(dat$CD3N_1), dat$VL_1, dat$A2_0, dat$A2_1)
  my1 = cbind(1, cd4draw1, cd8draw1, cd3draw1, VLdraw1, 1, 1) %*% as.numeric( Ys[(row+draw), 9:15] ) # always treated
  Ydraw1 = rnorm(n, mean = my1, sd = sqrt(Ys[(row+draw), 16]) )
  
  my0 = cbind(1, cd4draw0, cd8draw0, cd3draw0, VLdraw0, 0, 0) %*% as.numeric(Ys[(row+draw), 9:15]) # always treated
  Ydraw0 = rnorm(n, mean = my0, sd = sqrt(Ys[(row+draw), 16]) )
  
  draws11[,draw] = Ydraw1
  draws00[,draw] = Ydraw0
  
  dat$y11 = ifelse(dat$t1 == 1, log(dat$CD4N_2), Ydraw1)
  dat$y00 = ifelse(dat$t0 == 1, log(dat$CD4N_2), Ydraw0)

  diffs[draw] = sum(dat$y11 - dat$y00)/n
}

### pp for y11
draws11r = draws11[which(dat$t1 ==1),]
draws11d = data.frame(draws11r) %>% pivot_longer(cols = everything())

data = data.frame(name = "y", value = log(dat$CD4N_2)[which(dat$t1 ==1)])
plotdf1 = rbind(draws11d, data)

p1 = ggplot(plotdf1, aes(x=value, group = name) ) + geom_density(aes(color = name)) + 
  scale_color_manual(values=c( rep("cornflowerblue", ndraws), "firebrick4")) + 
  theme_minimal() + theme(legend.position = "none") + xlim(c(0,8.5))


### pp for y00
draws00r = draws00[which(dat$t0 ==1),]
draws00d = data.frame(draws00r) %>% pivot_longer(cols = everything())

data = data.frame(name = "y", value = log(dat$CD4N_2)[which(dat$t0 ==1)])
plotdf0 = rbind(draws00d, data)

p2 = ggplot(plotdf0, aes(x=value, group = name) ) + geom_density(aes(color = name)) + 
  scale_color_manual(values=c( rep("cornflowerblue", ndraws), "firebrick4")) + 
  theme_minimal() + theme(legend.position = "none") + xlim(c(0,8.5))

library(gridExtra)

g = arrangeGrob(p1, p2, nrow=1) #generates g
ggsave(file="whatever.pdf", g) #saves g





dat11 = filter(dat, t1 == 1) %>% dplyr::select(CD4N_2)
dat11$CD4N_2 = log(dat11$CD4N_2)
dat11.l = dat11 %>% pivot_longer(cols = c(y11, CD4N_2) )


ggplot(dat11.l, aes(x=value, group = name)) + 
  geom_density(aes(color = name))

hist(log(dat11$CD4N_2))
hist(dat11$y11)

library(latex2exp)
diffs = data.frame(diffs)
ggplot(diffs, aes(x=diffs)) + geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  geom_density(alpha=0.6) + theme_minimal(base_size = 25) + xlab(TeX('$\\phi$')) + 
  ggtitle(TeX('Posterior distribution of $\\phi$'))

ggsave("phipp.png", dpi = 320, width = 14, height = 5)

ggplot(diffs, aes(x=weight, color=sex, fill=sex)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  geom_density(alpha=0.6)

#par(mfrow = c(2,1))

#hist(Ydraw0, xlim = c(0,15))
#hist(Ydraw1, xlim = c(0,15))


# one pp check
pp = dat %>% select(CASEID, CD4N_2, t1, t0)
pp$CD4N_2 = log(pp$CD4N_2)
pp$y11 = Ydraw1
pp$y00 = Ydraw0
pp2 = filter(pp, t1==0) %>% select(CASEID, CD4N_2, y00) 
pp2l = pp2 %>% pivot_longer(!CASEID)

ggplot(pp2l, aes(value, fill = name)) + geom_density(alpha=0.3)

#qqplot(pp2$CD4N_2, pp2$y11)


#library(rstanarm)
#m0 = stan_glm(log(CD4N_2) ~ CD4N_1)
