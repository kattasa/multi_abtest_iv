
# Simulation study for causal using regularized I.V

# confounder parameters 
mu_u = 1
sd_u = 5

sd_x = 1
sd_y = 1

# cluster 
K = 50
n_k = 100
mu_z = 1  # randomized Norm design 
sd_z = 1:K # cluster-level variation 


# model parameters 
n = K*n_k
mu = 2
psi = 10
gamma = 10
beta = 1

# generate confounders 
z = rnorm(n_k,mu_z,sd_z[1])
indicator = rep(1,n_k)
for(k in 2:K){
  z_k = rnorm(n_k,mu_z,sd_z[k])
  z = c(z,z_k)
  indicator = c(indicator,rep(k,n_k))
}

u = rnorm(n,mean=mu_u,sd=sd_u) # un-obs confounder
x = z*mu+u*psi+rnorm(n,0,sd=sd_x) # obs confounder
y = beta*x + gamma*u + rnorm(n,0,sd=sd_y)

data = cbind(indicator,z,u,x,y)
data = data.frame(data)

# OLS 
OLS_model = lm(y~x,data=data)
OLS_beta = coef(OLS_model)[2]

# Two-stage least squares (TSLS)
x_fit = fitted(lm(x~z,data=data)) # stage 1 
TSLS_model = lm(y~x_fit,data=data) # stage 2 
TSLS_beta = coef(TSLS_model)[2] # second estimator is for beta 

# IVCV OLS 
IVCV_betas = numeric()
IVCV_loss = numeric()
IVCV_mse = numeric()

# IVCV TSLS 
TSLS_betas = numeric()
IVCV_TSLS_loss = numeric()
IVCV_TSLS_mse = numeric()

for(k in 1:K){
  # k = 1
  idx_k = data$indicator == k
  data_k = data[idx_k ,]
  
  # OLS 
  IVCV_model_k = lm(y~x,data=data_k)
  IVCV_betas[k] = coef(IVCV_model_k)[2]
  
  y_k_fit = fitted(IVCV_model_k)
  y_k = data_k$y
  x_k = data_k$x
  
  # IVCV_mse[k] = sum((y_k_fit-y_k)^2)
  IVCV_mse[k] = mean((y_k_fit - y_k)^2)
  IVCV_loss[k] = (mean(y_k) - mean(x)*IVCV_betas[k])^2
  
  # TSLS
  x_fit_k = fitted(lm(x~z,data=data_k))
  TSLS_model_k = lm(y~x_fit_k,data=data_k)
  TSLS_betas[k] = coef(TSLS_model_k)[2] # second estimator is for beta 
  y_k_fit = fitted(TSLS_model_k)
    
  IVCV_TSLS_loss[k] = (mean(y_k) - mean(x)*TSLS_betas[k])^2
  IVCV_TSLS_mse[k] = sum((y_k_fit-y_k)^2)
}

IVCV_Causal_Loss = ((IVCV_betas - beta))^2
IVCV_TSLS_Causal_Loss = ((TSLS_betas - beta))^2

# Normalized 
Normalization = function(x){ (x-min(x))/(max(x)-min(x))}

normalized_IVCV_mse = Normalization(x = IVCV_mse)
normalized_IVCV_loss = Normalization(x = IVCV_loss)
normalized_Causal_Loss = Normalization(x = Causal_Loss)

normalized_IVCV_TSLS_Causal_Loss = Normalization(x = IVCV_TSLS_Causal_Loss)
normalized_IVCV_TSLS_mse = Normalization(x = IVCV_TSLS_mse)
normalized_IVCV_TSLS_loss = Normalization(x = IVCV_TSLS_loss)

# estimation plot 
y_min = min(c(IVCV_betas,TSLS_betas))
y_max = max(c(IVCV_betas,TSLS_betas))
plot(sd_z,IVCV_betas,ylim=c(y_min,y_max),col=4,lwd=1,ylab="beta_hat") # IVCV_OLS
lines(sd_z,TSLS_betas,type="p",lwd=1) # IVCV_TSLS
abline(h=beta,col=2,lwd=3)
abline(h=OLS_beta,col=4,lwd=2,lty=2)
abline(h=TSLS_beta,col=1,lwd=2,lty=2)

# Loss plot 
plot(sd_z,normalized_IVCV_mse,col=4,lwd=1)
lines(sd_z,normalized_Causal_Loss,col=4,lwd=1)

plot(sd_z,normalized_IVCV_TSLS_mse,col=1,lwd=1,type="p")
lines(sd_z,normalized_IVCV_TSLS_Causal_Loss,col=1,lwd=1)
























































