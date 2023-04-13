#################################################################################
# Occupancy Power Analysis
# To recreate 3D plot, Figure 3
# Input: species-specific p and psi estimates, rest is simulated
# Output: 3D plot
# Notes: using code from Guillera-Arroita 2012, functions but modified
# Analysis and code from Robin Steenweg (Steenweg et al. In Review)
##################################################################################

	# ----------------------
	#  Load packages
	install.packages("rgl")
	library(plyr)
	library(rgl)

	
	
##################################################################################
#  Things to remember...
##################################################################################	

#  alpha: Type I error: : prob. of falsely detecting an increase or decrease in occupancy 
#   (i.e., false positive)

#  beta: Type II error: prob. of missing an increase or decrease in occupancy (i.e., false 
#   negative)

#  Power: 1- beta: prob. of correctly rejecting null hypothesis (null hypothesis is that 
#   occupancy is not changing).

#  Power is impacted by:
#   Effect size (i.e., magnitude of change in occupancy): value-based management decision 
#   Sampling size (how many surveys and sites)
#   Variance in occupancy and detection probability estimates (e.g, does the estimate of
#      occupancy probability vary from 0.2 to 0.9, or is it more narrow from 0.7 to 0.9?); 
#      larger variance makes changes harder to see; species-specific

#   psi: prob. of occupancy
#   p: prob. of detection
#   p*: cumulative prob. of detection (probability of detecting a species at least once across
#    all vists)
	
#  Word of wisdom from former grad student......start learning your probability 
#   distributions now:
# http://www.stat.umn.edu/geyer/old/5101/rlook.html#forward



##################################################################################
#  Define functions that calculate power (from Guillera-Arroita and Lahoz-Monfort 2012)
##################################################################################

	# ----------------------
	#  Function to simulate encounter history
	#   S: number of sites sampled
	#   K: number os sampling replicates (e.g., number of surveys)
	#   psi: prob. of occupancy
	#   p: prob. of detection
	
	simhist <- function(S, K, psi, p){
	
	# Set up a matrix to hold output from encounter history simulation loop
		h <- matrix(NA, nrow = S, ncol = K) 
		
		# Set up vector of "true" occupancy at each site
		z <- rep(0, n.sites) 
		# Ensures that the proportion of sites occupied is the same as specified psi.
		z[sample(1:n.sites, round(n.sites*psi), replace = F)] = 1 
		
		# Fill in empty matrix set up above
		for (ii in 1:S){
			#  Acutal encounter history is binomial draw ~ true occupancy and detection prob.
			h[ii,] <- rbinom(K, 1, p * z[ii]) 
			}
		
		# What we want given to use out of the equation
		return(h) 
		}

	# ----------------------
	#   Example runs of simhist()
	n.sites = 10
	ex1 <- simhist(S = 10, K = 3, psi = 0.6, p = 0.8)
	ex1
	ex2 <- simhist(S = 10, K = 3, psi = 0.6, p = 0.3)
	ex2
	
	# ----------------------
	#  Function that turns equation 3 from Steenweg paper into code
	#   S1 & S2: number of sites sampled for scenario 1 & scenario 2
	#   K1 & K2: number os sampling replicates (e.g., number of surveys) for scenario 1 & scenario 2
	#   psi1: prob. of occupancy for scenario 1 & scenario 2
	#   p1 & p2: prob. of detection for scenario 1 & scenario 2
	#   R: proportional decline in prob. of occupancy from scenario 1 to scenario 2
	#   alpha: Type I error level
	
	calcPowerFormula1tailed <- function(S1, S2, K1, K2, p1, p2, psi1, R, alpha){
		# Calculate psi2: determined by psi1 and R
		psi2 <- psi1 * (1 - R) 
		# Cumulative prob. of detection (p*) for p1 & p2 (equation 2 in paper)
		pp1 <- 1 - (1 - p1)^K1 
		pp2 <- 1 - (1 - p2)^K2 
		# Component of variance (second part of equation 1 from paper)
		F1 <- (1 - pp1)/(pp1 - K1 * p1 * (1 - p1)^(K1 - 1)) 
		F2 <- (1 - pp2)/(pp2 - K2 * p2 * (1 - p2)^(K2 - 1))
		#  Variance of prob. of occupancy estimates (full equation 1 from paper using
		#   second part calculated above - F1 & F2)
		var1 <- psi1 * (1 - psi1 + F1)/S1
		var2 <- psi2 * (1 - psi2 + F2)/S2
		
		# Equation 3 in paper...
		#   qnorm(): "quantile"; determine the Z-score of the pth quantile of the normal distribution? 
		#   For us, 5th percentile because we set alpha = 0.05
		limL <- (+qnorm(1 - alpha) * (sqrt(var1 + var2)) - (psi1 - psi2))/sqrt(var1 + var2)
		
		# Power...
		#   pnorm(): "probability", the cumulative distribution function (c. d. f.), inverse of qnorm()
		G <- 1 - pnorm(limL) # this is power
		
		# What we want given to use out of the equation
		return(G) 
		}

	# ----------------------
	#   Example runs of calcPowerFormula1tailed()
	ex3 <- calcPowerFormula1tailed(S1 = 65, S2 = 60, K1 = 4, K2 = 3, p1 =0.3, p2 = 0.5, psi1 = 0.7, R = 0.2, alpha = 0.05)

	ex3
	
	ex4 <- calcPowerFormula1tailed(S1 = 65, S2 = 60, K1 = 4, K2 =3, p1 = 0.8, p2 = 0.9, psi1 = 0.7, R = 0.2, alpha = 0.05)
	ex4
	
	ex5 <- calcPowerFormula1tailed(S1 = 65, S2 = 60, K1 = 4, K2 = 3, p1 = 0.3, p2 = 0.5, psi1 = 0.2, R = 0.2, alpha = 0.05)
	ex5
	

#######################################################################
#  Assess power across a variety of scenarios  (from Steenweg)
#######################################################################

	# ----------------------
	#  Function calculates power over a variety of scenarios
	# n.sites: number of sites sampled (S above)
	# n.rep.vect: vector of different number of power calculation replicates to calculate
	# psi.vect: vector of different occupancy probabilities (psi) to calculate
	# p.vect: vector of different detection probabilities (p) to calculate

	fnc.occu.pow <- function(n.sites, n.rep.vect, psi.vect, p.vect){
		# Starting iteration 
		xx = 0

		# Levels of prob. of occupancy decline (i.e., Effect size) to be assessed
		#decline.vect.max = c(0.05,0.1,0.15,0.2,0.25) 
		decline.vect.max = seq(0.01, 0.30, 0.001)

		# Total number of different scenarios of decline, psi and p that we will simulate
		total.sims = length(decline.vect.max) * length(psi.vect) * length(p.vect)
		
		# Matrix to hold the simulation outputs with 
		#  nrow = number of scenarios * number of power calculation replicates
		#  ncol = 4: (for eah of the set up values to be stored)
		stored = as.data.frame(matrix(NA, 
			nrow = length(n.rep.vect) * length(decline.vect.max) * length(psi.vect) * length(p.vect), 
			ncol = 4))
		#  Name columns
		names(stored) = c("psi", "p", "decline", "power")
		
		# Gives us a progress bar
		print(paste("Max simulations:", total.sims))
		
		# Use nested for loops to get every combination of the values of interest
		#  (i.e., the different scenarios) and replicate it the number of times desired (n.rep.vect)
		for(n.rep in n.rep.vect){ 
			for (psi in psi.vect){
				
				# Make sure the decline in prob. of occupancy can only be less than the 
				#  actual prob. of occupancy.
				decline.vect = decline.vect.max[1:length(which(decline.vect.max<(psi)))]
				
				for (p in p.vect){
					for (decline in decline.vect){
						
						# Move through loop from starting iteration (0) and making sure power 
						#  calculation was clear out from previous iteration
						xx=xx+1; pow=NULL
						
						#  Use function defined earlier with particular combination of values from 
						#  defined scenarios.
						pow = calcPowerFormula1tailed(S1 = n.sites, S2 = n.sites, K1 = n.rep, K2 = n.rep,
							p1 = p, p2 = p, psi1 = psi, R = decline/psi, alpha = 0.05)
						
						#  Store the outputs 
						stored$decline[xx] = decline
						stored$psi[xx] = psi
						stored$p[xx] = p
						stored$power[xx] = pow
						cat('\r', xx); flush.console()
						}
					}
				}
			}
		return(stored)
		}

	# ----------------------
	#   Example runs of fnc.occu.pow()
	ex6 <- fnc.occu.pow(n.sites = 60, n.rep.vect = c(5, 10, 15), psi.vect = c(0.3, 0.5, 0.7), 
		p.vect=c(0.1, 0.5, 0.7))
	dat_ex6 <- ex6 %>% mutate(group = paste(psi, p, sep = "_"))
	ex6_p <- ggplot() + geom_line(dat = dat_ex6, aes(x = decline, y = power, colour = group))
	ex6_p	

	
	
######################################################################## Large example run and 3D plot 
#######################################################################
#  Run power function 7 days; note S = n.sites, K = n.rep 
	dd <- fnc.occu.pow(n.sites = 183, n.rep.vect = 26, psi.vect = seq(0.1, 0.9, 0.05), p.vect = seq(0.05, 0.4, 0.05))

	#  Remove NA's at bottom because matrix too long
	dd = dd[!is.na(dd$decline),]

	#  Select rows that have the minimum decline but still have 80% power to detect
	#  decline = minimum value
	#  power >= .80
	#  Make a new data frame (power80), where you apply a function to select the minimum
	#   level of decline where power = 0.80 across each unique grouping of psi and p (i.e, 
	#   each scenario)
	power80 = ddply(dd, ~psi + p, function(df){
		ifelse(round(df[which.min(abs(df$power-0.80)),]$power,digits=2)>=0.8, 
			df[which.min(abs(df$power-0.80)),]$decline, 0)})
	names(power80)[3] = c("decline80")  
	save(power80, file = "power80.Rdata")
	
	
########################################################################  3D plot 
#######################################################################
	#  Format data for 3 dimenarions
	x = unique(power80$p)
	y = unique(power80$psi)
	z = t(as.matrix(tidyr::spread(power80, p, decline80)[,-1]))

	# ----------------------
	#  Make empty 3D plot with desired parameters
	plot3d(0, 0, 0, pch = "", xlab = "", ylab = "", zlab = "", axes = T, xlim = c(0, 0.4), ylim = c(0, 1), 
		zlim = c(0, 0.25), size = .5, type = "s")

	# ----------------------
	#  Set up colors based on z axis
	nbcol = 100
	color = rev(rainbow(nbcol, start = 0, end = 0.66))
	#color = gray.colors(nbcol, start = 0, end = 1) # for gray scale 
	zcol  = cut(z, nbcol)

	#  Add  the data to plot
	surface3d(unique(power80$p), unique(power80$psi), z, xlim = c(0, 0.4), ylim = c(0, 1), 
		zlim = c(0, 0.25), alpha = 0.8, col = color[zcol])

	#  Make background clean
	rgl.bbox(color = "grey50", # grey60 surface and black text
			 emission = "grey50", # emission color is grey50 
			 xlen = 0, ylen = 0, zlen = 0) # Don't add tick marks

	#  Set default color of future objects to black
	rgl.material(color = "black")
	
	#  Add axes to specific sides. Possible values are "x--", "x-+", "x+-", and "x++".
	axes3d(edge = "x--", ntick = 5, labels = c(0, 0.1, 0.2, 0.3, 1)) # cex to change font
	axes3d(edge = "y+-", ntick = 3, labels = c(0, 0.5, 1))
	axes3d(edge = "z--", ntick = 5, labels = c(0, "", 0.1, "", 0.2, ""))

	#  Add axis labels. 'line' specifies how far to set the label from the axis.
	rgl::mtext3d("Detection probability", edge = "x--", line = 2, at = 0.05)
	rgl::mtext3d("                        Occupancy probability", edge = "y+-", line = 2, at = 0.5)
	rgl::mtext3d("Absolute decline in", edge = "z--", line = 3.1, at = 0.155)
	rgl::mtext3d("occupancy detected", edge = "z--", line = 3, at = 0.14)
	rgl::mtext3d("with 80% power", edge = "z--", line = 3.3, at = 0.125)

	#  Add species-specific psi and p to 3D plot 
	multisp7d = data.frame(
		"species" = c('grizzly bear' ,'wolf' ,'mule deer' ,'moose' , 'white-tailed deer', 
			'black bear' ,'elk' ,'cougar' ,'lynx' ,'coyote' ,'red fox' ,'wolverine' ,'caribou'),
		"p" = c(0.196254497, 0.193955344, 0.181717663, 0.187934179, 0.244204191, 
			0.254281369, 0.237973916, 0.097674991, 0.115372418, 0.135632174, 
			0.154919324, 0.057619199, 0.17912606 ),
		"psi" = c(0.793857099, 0.683109834, 0.612291376, 0.603651072, 0.597613922, 
			0.579085479, 0.563764669, 0.291493677, 0.279054569, 0.255191376, 
			0.211412425, 0.19891626, 0.0351722),
		"decline80" = c(0.114, 0.126, 0.129, 0.129, 0.129, 0.129, 0.129, 0.116, 0.111,
			0.105, 0.096, 0.113, 0))
	spheres3d(multisp7d$p, multisp7d$psi, 0, radius = 0.012, fog = F)  
	
	#  Add lines
	for (j in c(1:13)){
	  segments3d(multisp7d$p[j], multisp7d$psi[j], c(0, multisp7d$decline80[j] + 0.01), lwd=2)
		}

	#  Add lines for named segments
	special = multisp7d[multisp7d$species %in% c("grizzly bear", "wolverine", "caribou"),]
	for (j in 1:3){ #
	  segments3d(special$p[j], special$psi[j], c(0, special$decline80[j] + 0.04), lwd=2)
	  text3d(special$p[j], special$psi[j], special$decline80[j]+0.05, special$species[j])
		}

		
		
########################################################################  Other power analysis options....
#######################################################################
#  R packages for simpler situations to calculate power analytically
install.packages("pwr")
#https://cran.r-project.org/web/packages/pwr/vignettes/pwr-vignette.html

#  Other examples to estimate power by simulating data
install.packages(c("paramtest", "nlme", "lavaan"))
# https://cran.r-project.org/web/packages/paramtest/vignettes/Simulating-Power.html

#  Examples in Bayesian framework
#  Introduction to WinBUGS for Ecologists: http://www.biometrica.tomsk.ru/lib/kery.pdf