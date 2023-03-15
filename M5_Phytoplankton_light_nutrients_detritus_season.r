########################### DAY 3 ###########################
library(deSolve)
library(fields)

# Create an empty list to add the parameters
param <- list()

# Define parameters 
param$I0 <- 1000 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
param$kw <- 0.0375 # light absorbed by water (1/m)
param$kc <- 0.05 # light absorbed by phyto plankton (m^2/mmol N)
param$gmax <- 1.5 # maximum growth (1/d)
param$alpha <- 0.1/(param$I0/200) # Slope of the PI-curve (1/(uEinstein*d)
param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
param$m <- 0.03 # specific loss rate (1/d)
param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
param$eps <- 0.1 # remineralization (1/d)
param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
param$ud <- 15 # Settling velocity (m/d)
param$up <- 0.5 # Settling velocity (m/d)
param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
param$dz <- 2.5 # grid spacing (m)
param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
param$n <- length(param$z) # number of grid cells
tee <- 5


# Define initial conditions - concentrations of phytoplankton and nutrients
PND <- c(rep(10,param$n),rep(param$Nb, param$n),rep(0,param$n))

# Create a function
# the 22*pi/365 to push it back to the 21st December
# lowest light at 21st of December 

CalLight <- function(t,P,D, param) {
  season <- 1+cos(pi + (2*pi)/365*t + 22*pi/365)
  
  Q <- param$kc * param$dz * ((cumsum(P) - P/2) + (cumsum(D) - D/2))
  
  I <- param$I0 * exp(-param$kw * param$z - Q) * season
  
  return(I)
}




# Create function that creates the differential equation at each step
FuncPND <- function(t, PND, param) {
  P <- PND[1:param$n]
  N <- PND[(1+param$n):(2*param$n)]
  D <- PND[(2*param$n+1):(3*param$n)]
 
  # Phytoplankton flux
  Jap <- rep(0,param$n+1)
  Jdp <- rep(0,param$n+1)
  
  for (i in 2:param$n){
    Jap[i] <- param$up * P[i-1]
    Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
  }
  
  # Advective flux boundary
  Jap[1] = 0
  Jap[param$n+1] = 0
  
  # Diffusive flux boundary
  Jdp[1] = 0
  Jdp[param$n+1] = 0
  
  Jp = Jap + Jdp
  
  # Nutrient flux
  Jdn <- rep(0,param$n+1)
  
  for (i in 2:param$n){
    Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
  }
  
  
  # Diffusive flux boundary
  Jdn[1] = 0
  Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
  
  Jn = Jdn
  
  # Detritus flux 
  Jad <- rep(0,param$n+1)
  Jdd <- rep(0,param$n+1)
  
  for (i in 2:param$n){
    Jad[i] <- param$ud * D[i-1]
    Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
  }
  
  # Advective flux boundary
  Jad[1] = 0
  Jad[param$n+1] = param$ud*D[param$n]
  
  # Diffusive flux boundary
  Jdd[1] = 0
  Jdd[param$n+1] = 0
  
  Jd = Jad + Jdd
  
  
  I <- CalLight(t,P,D,param)
  
  g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
  
  dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
  
  dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
  
  dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
  
  
  return(list(c(dPdT,dNdT, dDdT)))
}




# Define time step
time <-seq(0,365*3, by=1)

# Solve our differential equation
start.time <- Sys.time()
res <- ode(PND, time, FuncPND, param)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste("Runtime=",as.character(round(time.taken,1)),"s"))
P <- res[,2:(param$n+1)]
N <- res[,(param$n+2):(param$n*2+1)]
D <- res[,(param$n*2+2):(param$n*3+1)]



# Surface plot of phytoplankton
par(mfrow=c(1,3))
#par(mfrow=c(1,1))
image.plot(time/365, param$z, P, ylim = rev(range(param$z)), col = topo.colors(100)[20:100],
           xlab = "Time [days after January 1st]", ylab = "Depth [m]",
           main = "Phytoplankton [mmol N/m^3]")
box()


# Surface plot of Nutrients
image.plot(time/365, param$z, N, ylim = rev(range(param$z)), col = rev(rainbow(300, start= 0.1, end = 0.6)),
           xlab = "Time [days after January 1st]", ylab = "Depth [m]",
           main = "Nutrient [mmol N/m^3]")
box()


# Surface plot of Detritus
image.plot(time/365, param$z, D, ylim = rev(range(param$z)), col =rev(rainbow(300, start= 0.1, end = 0.6)),
           xlab = "Time [days after January 1st]", ylab = "Depth [m]",
           main = "Ditritus [mmol N/m^3]")
box()


par(mfrow=c(1,1))
plot(P[365+180,], param$z, ylim = rev(range(param$z)), pch = 19, col = "darkolivegreen", main = "1 days", ylab = "Depth", xlab = "Concentration")
points(N[365+180,], param$z, ylim = rev(range(param$z)), type = 'l', pch = 19, col = "red",cex = 0.5,lwd = 2)
points(D[365+180,], param$z, ylim = rev(range(param$z)), type = 'l', pch = 19, col = "blue",cex = 0.5,lwd = 2)


# VERTICAL PROFILE FOR STEADY STATE QUESTION - WINTER
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 0.3) 
plot(P[2*365,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkolivegreen", xlab = "", ylab = "Depth (m)", xaxt = "n")
points(D[2*365,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "blue", xlab = "", ylab = "")
axis(1, col = "black", col.ticks = "black", col.axis = "black", col.lab="black")
mtext("Phytoplankton and detritus concentration [mmolN/m^3]", side = 1, line = 3, col = "black") 
par(new = TRUE) 
plot(N[2*365,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkred", 
     axes = FALSE, xlab = "", ylab = "")
axis(side = 3, at = pretty(range(N[2*365,])), col = "darkred", col.ticks = "darkred", col.axis = "darkred", col.lab="darkred")
mtext("Nutrient concentration [mmolN/m^3]", side = 3, line = 3, col = "darkred") 
legend(18,4, legend = c("Phytoplankton","Detritus","Nutrients"), col = c("darkolivegreen","blue","darkred"), lty = 1, cex = 0.9, lwd=3, box.col = "white")
legend(8,5, legend = substitute(paste(bold("1st of January"))), box.col = "white")


# VERTICAL PROFILE FOR STEADY STATE QUESTION - SPRING
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 0.3) 
plot(P[2*365+90,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkolivegreen", xlab = "", ylab = "Depth (m)", xaxt = "n")
points(D[2*365+90,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "blue", xlab = "", ylab = "")
axis(1, col = "black", col.ticks = "black", col.axis = "black", col.lab="black")
mtext("Phytoplankton and detritus concentration [mmolN/m^3]", side = 1, line = 3, col = "black") 
par(new = TRUE) 
plot(N[2*365+90,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkred", 
     axes = FALSE, xlab = "", ylab = "")
axis(side = 3, at = pretty(range(N[2*365+90,])), col = "darkred", col.ticks = "darkred", col.axis = "darkred", col.lab="darkred")
mtext("Nutrient concentration [mmolN/m^3]", side = 3, line = 3, col = "darkred") 
#legend("right", legend = c("Phytoplankton","Detritus","Nutrients"), col = c("darkolivegreen","blue","darkred"), lty = 1, cex = 0.8, lwd=3)
legend(8,5, legend = substitute(paste(bold("1st of April"))),box.col = "white")

# VERTICAL PROFILE FOR STEADY STATE QUESTION - SUMMER
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 0.3) 
plot(P[2*365+181,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkolivegreen", xlab = "", ylab = "Depth (m)", xaxt = "n")
points(D[2*365+181,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "blue", xlab = "", ylab = "")
axis(1, col = "black", col.ticks = "black", col.axis = "black", col.lab="black")
mtext("Phytoplankton and detritus concentration [mmolN/m^3]", side = 1, line = 3, col = "black") 
par(new = TRUE) 
plot(N[2*365+181,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkred", 
     axes = FALSE, xlab = "", ylab = "")
axis(side = 3, at = pretty(range(N[2*365+181,])), col = "darkred", col.ticks = "darkred", col.axis = "darkred", col.lab="darkred")
mtext("Nutrient concentration [mmolN/m^3]", side = 3, line = 3, col = "darkred") 
#legend("right", legend = c("Phytoplankton","Detritus","Nutrients"), col = c("darkolivegreen","blue","darkred"), lty = 1, cex = 0.8, lwd=3)
legend(8,5, legend = substitute(paste(bold("1st of July"))),  box.col = "white")


# VERTICAL PROFILE FOR STEADY STATE QUESTION - FALL
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 0.3) 
plot(P[2*365+273,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkolivegreen", xlab = "", ylab = "Depth (m)", xaxt = "n")
points(D[2*365+273,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "blue", xlab = "", ylab = "")
axis(1, col = "black", col.ticks = "black", col.axis = "black", col.lab="black")
mtext("Phytoplankton and detritus concentration [mmolN/m^3]", side = 1, line = 3, col = "black") 
par(new = TRUE) 
plot(N[2*365+273,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkred", 
     axes = FALSE, xlab = "", ylab = "")
axis(side = 3, at = pretty(range(N[2*365+273,])), col = "darkred", col.ticks = "darkred", col.axis = "darkred", col.lab="darkred")
mtext("Nutrient concentration [mmolN/m^3]", side = 3, line = 3, col = "darkred") 
#legend("right", legend = c("Phytoplankton","Detritus","Nutrients"), col = c("darkolivegreen","blue","darkred"), lty = 1, cex = 0.8, lwd=3)
legend(8,5, legend = substitute(paste(bold("1st of October"))), box.col = "white")




# LIMITATION - WINTER
Ieq <- CalLight(2*365,P[2*365,],D[2*365,], param)
fi <- param$alpha*Ieq / sqrt((param$alpha*Ieq)**2 + param$gmax**2)
fn <- N[2*365,] / (N[2*365,] + param$HN)
par(mar = c(5, 4, 4, 4) + 0.3) 
plot(fi, param$z, ylim = rev(range(param$z)), xlim=c(0,1),type='l', lwd=3, col = "gold", xlab = "Unit interval (Unitless)", ylab = "Depth (m)", main="Limiting factor for phytoplankton growth")
points(fn, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "brown")
legend (0.7,30, legend = c("Light", "Nutrients"),lwd=3,col=c("gold","brown"), cex = 1.1, box.col = "white")
legend(0.3,4, legend = substitute(paste(bold("1st of January"))),  box.col = "white")

# LIMITATION - SPRING
Ieq <- CalLight(2*365+90,P[2*365+90,],D[2*365+90,], param)
fi <- param$alpha*Ieq / sqrt((param$alpha*Ieq)**2 + param$gmax**2)
fn <- N[2*365+90,] / (N[2*365+90,] + param$HN)
par(mar = c(5, 4, 4, 4) + 0.3) 
plot(fi, param$z, ylim = rev(range(param$z)), xlim=c(0,1),type='l', lwd=3, col = "gold", xlab = "Unit interval (Unitless)", ylab = "Depth (m)", main="Limiting factor for phytoplankton growth")
points(fn, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "brown")
#legend (0.7,30, legend = c("Light", "Nutrients"),lwd=3,col=c("gold","brown"), cex = 1.0,box.col = "white")
legend(0.3,4, legend = substitute(paste(bold("1st of April"))),  box.col = "white")

# LIMITATION - SUMMER
Ieq <- CalLight(2*365+181,P[2*365+181,],D[2*365+181,], param)
fi <- param$alpha*Ieq / sqrt((param$alpha*Ieq)**2 + param$gmax**2)
fn <- N[2*365+181,] / (N[2*365+181,] + param$HN)
par(mar = c(5, 4, 4, 4) + 0.3) 
plot(fi, param$z, ylim = rev(range(param$z)), xlim=c(0,1),type='l', lwd=3, col = "gold", xlab = "Unit interval (Unitless)", ylab = "Depth (m)", main="Limiting factor for phytoplankton growth")
points(fn, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "brown")
#legend (0.7,30, legend = c("Light", "Nutrients"),lwd=3,col=c("gold","brown"), cex = 1.0,box.col = "white")
legend(0.3,4, legend = substitute(paste(bold("1st of July"))),  box.col = "white")

# LIMITATION - FALL
Ieq <- CalLight(2*365+273,P[2*365+273,],D[2*365+273,], param)
fi <- param$alpha*Ieq / sqrt((param$alpha*Ieq)**2 + param$gmax**2)
fn <- N[2*365+273,] / (N[2*365+273,] + param$HN)
par(mar = c(5, 4, 4, 4) + 0.3) 
plot(fi, param$z, ylim = rev(range(param$z)), xlim=c(0,1),type='l', lwd=3, col = "gold", xlab = "Unit interval (Unitless)", ylab = "Depth (m)", main="Limiting factor for phytoplankton growth")
points(fn, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "brown")
#legend (0.7,30, legend = c("Light", "Nutrients"),lwd=3,col=c("gold","brown"), cex = 1.0,box.col = "white")
legend(0.3,4, legend = substitute(paste(bold("1st of October"))),  box.col = "white")
