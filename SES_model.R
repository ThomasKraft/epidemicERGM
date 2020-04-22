######## SES ##########

# load packages
library(statnet)
library(EpiModel)
library(ndtv)

######## SES ##########
# 1) Network with SES
# initialize empty network with desired attributes
net <- network.initialize(n = 200, directed = FALSE)
net <- set.vertex.attribute(net, "SES", rep(0:1, each = 100))

# Setting model parameters
# NOTE: this is where we can add attributes that will define structure in the way we like based on an ERGM formula
formation <- ~edges + nodematch("SES", diff=TRUE) + nodefactor("SES") 
target.stats <- c(150, 30, 15, 20) # correspond to formation formula above

# because networks generated are dynamic, we also set parameters to determine the rate of dissolution of ties
coef.diss <- dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("SES", diff=TRUE)), 
                               duration = c(1, 10, 10))

# Fitting model with desired parameters
mod1 <- netest(net, formation, target.stats, coef.diss, edapprox = T)

# Model diagnostics
dx <- netdx(mod1, nsims = 5, nsteps = 200, dynamic=F)
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(dx)


# Simulate epidemic
param <- param.net(inf.prob = 0.3, act.rate = 1, rec.rate=.05)
status.vector <- rbinom(200, 1, 0.03)   # seed initial infected individuals
status.vector <- ifelse(status.vector == 1, "i", "s")
init <- init.net(status.vector = status.vector)

# Structural model features
control <- control.net(type = "SIR", nsteps = 300, nsims = 10, epi.by = "SES", ncores=2)

# simulate!
sim1 <- netsim(mod1, param, init, control)

#### WHAT CAN WE EXTRACT FROM THE "sim1" OBJECT ?? ####




# # The rest is just plotting
# # Plots of how the epidemic spread
plot(sim1, mean.line = FALSE, qnts = FALSE, sim.lines = TRUE)
# 
# # Plots of how the simulated dynamic networks looked at different timepoints
# par(mfrow = c(1,2), mar = c(0,0,1,0))
# plot(sim1, type = "network", at = 1, col.status = TRUE,
#      main = "Prevalence at t1")
# plot(sim1, type = "network", at = 50, col.status = TRUE,
#      main = "Prevalence at t100")
# 
# 
# # plot networks at specific points showing race
nw <- get_network(sim1, sim = 1)
out <- network.extract(nw, at=20)
plot(out, vertex.col="SES")

# Make an animated plot of the networks over a specific duration
slice.par<-list(start=1, end=30, interval=1, aggregate.dur=1,rule="latest")
compute.animation(nw,slice.par=slice.par)

render.d3movie(nw,
               filename="SES_network.html",
               edge.col="darkgray",displaylabels=TRUE,
               label.cex=.6,label.col="blue",
               output.mode = 'HTML',
               vertex.col="SES")

render.animation(nw, vertex.cex=0.6, vertex.col="SES")
saveVideo(ani.replay(),video.name="SES_network.mp4",
          other.opts="-b 5000k",clean=TRUE)



######## RACE ##########
# 1) Network with RACE
# initialize empty network with desired attributes
net <- network.initialize(n = 200, directed = FALSE)
net <- set.vertex.attribute(net, "race", c(rep(0:1, each = 70), rep(2, 60)))

# Setting model parameters
# NOTE: this is where we can add attributes that will define structure in the way we like based on an ERGM formula
formation <- ~edges + nodematch("race")
target.stats <- c(200, 180) # correspond to formation formula above

# stats_convert <- function(n,mean.degree,n_race_1,mean_deg_race_1,prcnt_homophil) {
#   return(c(mean.degree*n/2,mean_deg_race_1*n_race_1/2, prcnt_homophil*mean.degree*n/200,n/10 ))  
# }
# target.stats <- stats_convert (450, 1, 150, 1, 100)

# because networks generated are dynamic, we also set parameters to determine the rate of dissolution of ties
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 5)

# Fitting model with desired parameters
mod1 <- netest(net, formation, target.stats, coef.diss, edapprox = T)


# Simulate epidemic
param <- param.net(inf.prob = 0.3, act.rate = 1, rec.rate=0.001)
status.vector <- rbinom(200, 1, 0.03)   # seed initial infected individuals
status.vector <- ifelse(status.vector == 1, "i", "s")
init <- init.net(status.vector = status.vector)

# Structural model features
control <- control.net(type = "SIR", nsteps = 200, nsims = 10, epi.by = "race", ncores=2)

# simulate!
sim1 <- netsim(mod1, param, init, control)


# Plots of how the epidemic spread
plot(sim1, mean.line = FALSE, qnts = FALSE, sim.lines = TRUE)

# Plots of how the simulated dynamic networks looked at different timepoints
par(mfrow = c(1,2), mar = c(0,0,1,0))
plot(sim1, type = "network", at = 1, col.status = TRUE,
     main = "Prevalence at t1")
plot(sim1, type = "network", at = 50, col.status = TRUE,
     main = "Prevalence at t100")


# plot networks at specific points showing race
nw <- get_network(sim1, sim = 1)
out <- network.extract(nw, at=100)
plot(out, vertex.col="race")

# Make an animated plot of the networks over a specific duration
slice.par<-list(start=1, end=10, interval=1, aggregate.dur=1,rule="latest")
compute.animation(nw,slice.par=slice.par)

render.d3movie  (nw,
                 filename="homophily_by_race_network.html",
                 edge.col="darkgray",displaylabels=TRUE,
                 label.cex=.6,label.col="blue",
                 output.mode = 'inline',
                 )

render.animation(nw, vertex.col="race")
saveVideo(ani.replay(),video.name="homophily_by_race_network.mp4",
          other.opts="-b 5000k",clean=TRUE)
