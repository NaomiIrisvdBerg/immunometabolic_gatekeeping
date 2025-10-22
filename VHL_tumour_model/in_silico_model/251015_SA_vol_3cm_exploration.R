##The code below contains the mathematical model presented in Supplementary Materials III: Mathematical model exploring role of aerobic glycolysis in 3cm-rule
#it contains the code needed to run the model and generate associated figures as presented in Suppl. Fig. 7.
##note that in this code, we refer to the Immune Dysfunction Score (IDS) as the Metastatic Potential Index (MPI). This was renamed in the manuscript but is quantified the same way.

##force R to install only binary versions:
#options(pkgType = "binary")

##Install from CRAN (with binaries only if uncommented above)
install.packages(c("munsell", "deSolve", "ggplot2", "dplyr", "tidyr", "stringi", "purrr", "reshape2", "tidyverse", "car", "ggpubr", "vegan"))

##load the libraries
library(munsell) ##we ran the code with v 0.5.1
library(deSolve) ##v 1.40
library(ggplot2) ##v 4.0.0
library(dplyr) ##v 1.1.4
library(tidyr) ##v 1.3.1
library(reshape2) ##v 1.4.4
library(purrr) ##v 1.1.0
library(car) ##v 3.1.3
library(ggpubr) ##v 0.6.1

install.packages("ecoCopula")
library(ecoCopula) ##v 1.0.2
install.packages("dismo")
library(dismo) ##v 1.3.16

install.packages("spatstat")
library(spatstat) ##v 3.4.1


##Note that the below analysis explores a hypothesis that relies on the assumption that the metastatic potential of VHL renal tumours is negated by the tumour TME;
#in particular; we assume that the metastatic potential of a tumour is limited for as long as the tumour is in equilibrium with a functioning immune system
#since the metastatic potential of tumours is consistently observed at a given size threshold (i.e., the "3cm rule"), irrespective of a particular genetic profile of a tumour,
#we here assume that: 
#A. the metastatic potential is negated by a functioning immune system; and 
#B. that an exhausted immune system thus translates in an enhanced metastatic potential of the tumour --> thus an earlier inflection point in the tumour's metastatic potential as a function of size
#we explore that immune exhaustion is enhanced by accumulation of metabolic byproducts resulting from the tumour's (inefficient) metabolism/upregulated glycolysis
#which in turn we assume scales nonlinearly with tumour volume
#and thus that the inflection point observed in tumour metastatic potential as a function of tumour size (i.e., the "3cm rule") reflects an immuno-metabolic tipping point. 




####No suppression control:
###Sim func where we draw params from biologically informed/approximated distributions
simulate_profiles <- function(R, N = 400) {
  x <- seq(1e-4, R, length.out = N)
  dx <- diff(x)[1]
  
  #Basic geometry; how radius scales with surface area & vol.
  vol <- (4/3) * pi * R^3
  surf <- 4 * pi * R^2
  SA_V <- surf / vol  #surface area to volume ratio 
  
  #Params
  ##production/consumption per unit volume
  p_L0 <- 1     ##lactate produced per unit tumour volume
  p_T0 <- 1     ##O2 consumed per unit volume
  p_I0 <- 1     ##O2 consumed when immune cells are present at I = 1 (locally)
  
  ##boundary conditions scale with SA:Vol
  I_ext_base <- 1 ##functioning immune availability at tumour periphery, normalised to 1 to represent functioing immune presence at tumour edge
  O_ext_base <- 1 #normalised; assumed normoxic interface with surrounding env. (e.g, vasculature)
  L_ext_base <- pmax(0, 1 - SA_V) ##poorer lactate clearance with decreasing SA:Vol; captures the idea that larger tumours cannot clear lactate as efficiently, even at the boundary
  
  #diffusion const. (fixed)
  D_L <- 1 #normalised, baseline
  D_I <- 0.01 #immune cells infiltrate much slower/need active migration as they are bigger (constrained e.g., by stroma)
  D_O <- 10 #diffusion of O2 has to be at least one order of magnitude higher than those of larger species lactate or immune cells
  
  #Immune suppression by lactate ( Lactic acid inhibits T cell motility and cytokine release + induces transcriptional suppression of effector programs)
  #Lactate suppresses CTL activation; reversible by pH buffering
  a_L <- 0#runif(1, 0.1, 0.9) #Captures immunosuppressive effects of low pH / lactic acidosis on T cells, via e.g., acid-sensitive signaling or T cell arrest (10-90% of normal)
  a_O <- 0#runif(1, 0.1, 0.9) #captures T cell exhaustion/suppression from (real) hypoxia (10-90% of normal activity)
  
  ##Thresholds for MPI
  O_thresh <- runif(1, 0.1, 0.5) #hypoxia is often defined within the range of 10%-50% of O2 conc. of normoxia (i.e., normal oxygen conc in healthy tissue)
  I_thresh <- 0.2 #Immune cell density below which tumour is considered "immune suppressed" for our score (MPI)
  ##so immune exhaustion when we have only 0.2 of the functioning immune cell to tumour cell saturation/equilibrium state (i..e, 1)
  
  ##Scaled params
  p_L <- p_L0 * vol ##lactate production scales with total metabolic demand, which grows with number of cells = volume
  p_T <- p_T0 * vol ##O2 consumption scales with total metabolic demand, which grows with number of cells = volume
  p_I <- p_I0 #fixed rate. While for tumour, total consumption scales with volume, for immune cells: localised consumption scales with local density
  I_ext <- I_ext_base * SA_V #Immune cells enter tumour from the outside (e.g., via lymph or vasculature), so immune access is also constrained by surface area to volume ratio
  O_ext <- O_ext_base * SA_V #Oxygen enters tumour via surface (e.g., blood vessels or capillary exchange). As volume increases faster than surface area, the same surface must serve a larger mass
  
  ##initialise profiles
  L <- rep(0, N) #assumes no pre-acidification
  I <- rep(0.1, N) #represents gradual infiltration or boundary supply
  O <- rep(1, N) #high initial value needed to allow depletion toward core
  
  #Looped solve
  tol <- 1e-4
  for (iter in 1:1000) {
    L_old <- L
    I_old <- I
    O_old <- O
    
    for (i in 2:(N-1)) {
      r <- x[i]
      
      #Lactate diffusion + production (discretised steady-state spherical diffusion eq.)
      L[i] <- ((r + dx/2)^2 * L[i+1] + (r - dx/2)^2 * L[i-1] + dx^2 * p_L * r^2 / D_L) /
        ((r + dx/2)^2 + (r - dx/2)^2)
      #^ lactate diffuses from production sources (tumour cells) 
      #^uniform production throughout the tumour (irregardless of local O2 conc), scaled by tumour volume
      #^lactate then diffuses outward but clearance is harder in larger tumours due to higher volume but falling SA
      
      ##Immune diffusion with suppression
      hypoxia_penalty <- ifelse(O[i] < O_thresh, 1, 0)
      denom <- ((r + dx/2)^2 + (r - dx/2)^2 + dx^2 * r^2 / D_I * (a_L * L[i] + a_O * hypoxia_penalty))
      #denom <- ((r + dx/2)^2 + (r - dx/2)^2 + dx^2 * a_L * L[i] * r^2 / D_I) #if not accounting for immune cell exhaustion
      #denom <- ((r + dx/2)^2 + (r - dx/2)^2 + dx^2 * r^2 / D_I * (a_L * L[i] + a_O * (1 - O[i])))  #if accounting for exhaustion smoothly
      I[i] <- ((r + dx/2)^2 * I[i+1] + (r - dx/2)^2 * I[i-1]) / denom
      #^immune cells diffuse inward from the tumour boundary (where we assume functioning vasculature)
      #^movement suppressed by high lactate conc. (acidosis), simulating T cell exhaustion/dysfunction
      #^if a_L is 0, this reduces to no lactate-mediated immune exhaustion (implicit suppression for more numerical stability)
      #^when O2 falls below threshold (O_thresh), hypoxia may occur, which can exhaust immune cells
      #Note that I[x] is a dimensionless/normalised concentration representing relative immune cell density at location x
      #It is not discrete, as we define it as a fractional occupancy (where 1 = maximum functioning immune infiltration [representing e.g., 1 immune cell per tumour cell], 0 = no immune infiltration)
      
      
      
      # Oxygen diffusion with consumption (oxygen consumption from immune cells, proportional to their local presence, and oxygen consumption cancer cells scaling with total volume)
      total_consump <- p_T + p_I * I[i] #oxygen demand from immune cells is local, and scaled by relative density
      O[i] <- ((r + dx/2)^2 * O[i+1] + (r - dx/2)^2 * O[i-1] - dx^2 * total_consump * r^2 / D_O) /
        ((r + dx/2)^2 + (r - dx/2)^2)
      ##O2 is supplied from the outside (tumour boundary/surface area), as we assume functional vasculature there
      #it diffuses inward but is consumed by both tumour cells and infiltrating immune cells
    }
    
    ##boundary conditions
    L[N] <- L_ext_base  ##if 0, we assume lactate is rapidly cleared at the tumour edge regardless of size
    L[1] <- L[2]
    I[N] <- I_ext
    I[1] <- I[2]
    O[N] <- O_ext
    O[1] <- O[2]
    
    if (max(abs(L - L_old), abs(I - I_old), abs(O - O_old)) < tol) break
  }
  
  ##MPI zone: when functioning immune cells cannot infiltrate accordingly/ infiltrated cells become exhausted/dysfunctional
  MPI_zone <- which(I < I_thresh)
  MPI_frac <- length(MPI_zone) / N
  
  tibble(
    x = x, L = L, I = I, O = O,
    MPI = x %in% x[MPI_zone],
    radius = R,
    MPI_fraction = MPI_frac
  )
}

#sim settings
radii <- seq(0.5, 5, by = 0.25)
n_reps <- 50  #param draws per tested radius

#run
set.seed(42)
all_results <- map_dfr(radii, function(r) {
  map_dfr(1:n_reps, function(rep) {
    simulate_profiles(R = r) %>%
      mutate(sim_id = rep)
  })
})

##conc profiles
results_long <- all_results %>%
  tidyr::pivot_longer(cols = c(L, I, O), names_to = "Component", values_to = "Value")

ggplot(results_long, aes(x = x, y = Value, color = Component)) +
  geom_line(linewidth = 1) +
  facet_wrap(~radius, labeller = label_both, scales = "free_x") +
  labs(
    title = "Concentration profiles vs tumour radius",
    x = "Radial distance (cm)", y = "value"
  ) +
  theme_minimal()

##MPI emergence (not sure what else to call this)
mpi_summary <- all_results %>%
  dplyr::group_by(radius) %>%
  dplyr::summarise(
    median_MPI = median(MPI_fraction),
    lower = quantile(MPI_fraction, 0.25),
    upper = quantile(MPI_fraction, 0.75)
  )


ggplot(mpi_summary, aes(x = radius, y = median_MPI)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(
    name = "Tumour radius (cm)",
    sec.axis = sec_axis(~ . * 2, name = "Tumour diameter (cm)")
  ) +
  labs(
    title = "No suppression control",
    y = "Median MPI fraction"
  ) +
  ylim(0.7,1) +
  theme_minimal()










####T cell exhaustion/suppression from hypoxia only:
###Sim func where we draw params from biologically informed/approximated distributions
simulate_profiles <- function(R, N = 400) {
  x <- seq(1e-4, R, length.out = N)
  dx <- diff(x)[1]
  
  #Basic geometry; how radius scales with surface area & vol.
  vol <- (4/3) * pi * R^3
  surf <- 4 * pi * R^2
  SA_V <- surf / vol  #surface area to volume ratio 
  
  #Params
  ##production/consumption per unit volume
  p_L0 <- 1     ##lactate produced per unit tumour volume
  p_T0 <- 1     ##O2 consumed per unit volume
  p_I0 <- 1     ##O2 consumed when immune cells are present at I = 1 (locally)
  
  ##boundary conditions scale with SA:Vol
  I_ext_base <- 1 ##functioning immune availability at tumour periphery, normalised to 1 to represent functioing immune presence at tumour edge
  O_ext_base <- 1 #normalised; assumed normoxic interface with surrounding env. (e.g, vasculature)
  L_ext_base <- pmax(0, 1 - SA_V) ##poorer lactate clearance with decreasing SA:Vol; captures the idea that larger tumours cannot clear lactate as efficiently, even at the boundary
  
  #diffusion const. (fixed)
  D_L <- 1 #normalised, baseline
  D_I <- 0.01 #immune cells infiltrate much slower/need active migration as they are bigger (constrained e.g., by stroma)
  D_O <- 10 #diffusion of O2 has to be at least one order of magnitude higher than those of larger species lactate or immune cells
  
  #Immune suppression by lactate ( Lactic acid inhibits T cell motility and cytokine release + induces transcriptional suppression of effector programs)
  #Lactate suppresses CTL activation; reversible by pH buffering
  a_L <- 0#runif(1, 0.1, 0.9) #Captures immunosuppressive effects of low pH / lactic acidosis on T cells, via e.g., acid-sensitive signaling or T cell arrest (10-90% of normal)
  a_O <- runif(1, 0.1, 0.9) #captures T cell exhaustion/suppression from (real) hypoxia (10-90% of normal activity)
  
  ##Thresholds for MPI
  O_thresh <- runif(1, 0.1, 0.5) #hypoxia is often defined within the range of 10%-50% of O2 conc. of normoxia (i.e., normal oxygen conc in healthy tissue)
  I_thresh <- 0.2 #Immune cell density below which tumour is considered "immune suppressed" for our score (MPI)
  ##so immune exhaustion when we have only 0.2 of the functioning immune cell to tumour cell saturation/equilibrium state (i..e, 1)
  
  ##Scaled params
  p_L <- p_L0 * vol ##lactate production scales with total metabolic demand, which grows with number of cells = volume
  p_T <- p_T0 * vol ##O2 consumption scales with total metabolic demand, which grows with number of cells = volume
  p_I <- p_I0 #fixed rate. While for tumour, total consumption scales with volume, for immune cells: localised consumption scales with local density
  I_ext <- I_ext_base * SA_V #Immune cells enter tumour from the outside (e.g., via lymph or vasculature), so immune access is also constrained by surface area to volume ratio
  O_ext <- O_ext_base * SA_V #Oxygen enters tumour via surface (e.g., blood vessels or capillary exchange). As volume increases faster than surface area, the same surface must serve a larger mass
  
  ##initialise profiles
  L <- rep(0, N) #assumes no pre-acidification
  I <- rep(0.1, N) #represents gradual infiltration or boundary supply
  O <- rep(1, N) #high initial value needed to allow depletion toward core
  
  #Looped solve
  tol <- 1e-4
  for (iter in 1:1000) {
    L_old <- L
    I_old <- I
    O_old <- O
    
    for (i in 2:(N-1)) {
      r <- x[i]
      
      #Lactate diffusion + production (discretised steady-state spherical diffusion eq.)
      L[i] <- ((r + dx/2)^2 * L[i+1] + (r - dx/2)^2 * L[i-1] + dx^2 * p_L * r^2 / D_L) /
        ((r + dx/2)^2 + (r - dx/2)^2)
      #^ lactate diffuses from production sources (tumour cells) 
      #^uniform production throughout the tumour (irregardless of local O2 conc), scaled by tumour volume
      #^lactate then diffuses outward but clearance is harder in larger tumours due to higher volume but falling SA
      
      ##Immune diffusion with suppression
      hypoxia_penalty <- ifelse(O[i] < O_thresh, 1, 0)
      denom <- ((r + dx/2)^2 + (r - dx/2)^2 + dx^2 * r^2 / D_I * (a_L * L[i] + a_O * hypoxia_penalty))
      #denom <- ((r + dx/2)^2 + (r - dx/2)^2 + dx^2 * a_L * L[i] * r^2 / D_I) #if not accounting for immune cell exhaustion
      #denom <- ((r + dx/2)^2 + (r - dx/2)^2 + dx^2 * r^2 / D_I * (a_L * L[i] + a_O * (1 - O[i])))  #if accounting for exhaustion smoothly
      I[i] <- ((r + dx/2)^2 * I[i+1] + (r - dx/2)^2 * I[i-1]) / denom
      #^immune cells diffuse inward from the tumour boundary (where we assume functioning vasculature)
      #^movement suppressed by high lactate conc. (acidosis), simulating T cell exhaustion/dysfunction
      #^if a_L is 0, this reduces to no lactate-mediated immune exhaustion (implicit suppression for more numerical stability)
      #^when O2 falls below threshold (O_thresh), hypoxia may occur, which can exhaust immune cells
      #Note that I[x] is a dimensionless/normalised concentration representing relative immune cell density at location x
      #It is not discrete, as we define it as a fractional occupancy (where 1 = maximum functioning immune infiltration [representing e.g., 1 immune cell per tumour cell], 0 = no immune infiltration)
      
      
      
      # Oxygen diffusion with consumption (oxygen consumption from immune cells, proportional to their local presence, and oxygen consumption cancer cells scaling with total volume)
      total_consump <- p_T + p_I * I[i] #oxygen demand from immune cells is local, and scaled by relative density
      O[i] <- ((r + dx/2)^2 * O[i+1] + (r - dx/2)^2 * O[i-1] - dx^2 * total_consump * r^2 / D_O) /
        ((r + dx/2)^2 + (r - dx/2)^2)
      ##O2 is supplied from the outside (tumour boundary/surface area), as we assume functional vasculature there
      #it diffuses inward but is consumed by both tumour cells and infiltrating immune cells
    }
    
    ##boundary conditions
    L[N] <- L_ext_base  ##if 0, we assume lactate is rapidly cleared at the tumour edge regardless of size
    L[1] <- L[2]
    I[N] <- I_ext
    I[1] <- I[2]
    O[N] <- O_ext
    O[1] <- O[2]
    
    if (max(abs(L - L_old), abs(I - I_old), abs(O - O_old)) < tol) break
  }
  
  ##MPI zone: when functioning immune cells cannot infiltrate accordingly/ infiltrated cells become exhausted/dysfunctional
  MPI_zone <- which(I < I_thresh)
  MPI_frac <- length(MPI_zone) / N
  
  tibble(
    x = x, L = L, I = I, O = O,
    MPI = x %in% x[MPI_zone],
    radius = R,
    MPI_fraction = MPI_frac
  )
}

#sim settings
radii <- seq(0.5, 5, by = 0.25)
n_reps <- 50  #param draws per tested radius

#run
set.seed(42)
all_results <- map_dfr(radii, function(r) {
  map_dfr(1:n_reps, function(rep) {
    simulate_profiles(R = r) %>%
      mutate(sim_id = rep)
  })
})

##conc profiles
results_long <- all_results %>%
  pivot_longer(cols = c(L, I, O), names_to = "Component", values_to = "Value")

ggplot(results_long, aes(x = x, y = Value, color = Component)) +
  geom_line(linewidth = 1) +
  facet_wrap(~radius, labeller = label_both, scales = "free_x") +
  labs(
    title = "Concentration profiles vs tumour radius",
    x = "Radial distance (cm)", y = "value"
  ) +
  theme_minimal()

##MPI emergence (not sure what else to call this)
mpi_summary <- all_results %>%
  group_by(radius) %>%
  summarise(
    median_MPI = median(MPI_fraction),
    lower = quantile(MPI_fraction, 0.25),
    upper = quantile(MPI_fraction, 0.75)
  )


ggplot(mpi_summary, aes(x = radius, y = median_MPI)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(
    name = "Tumour radius (cm)",
    sec.axis = sec_axis(~ . * 2, name = "Tumour diameter (cm)")
  ) +
  labs(
    title = "T cell exhaustion from hypoxia",
    y = "Median MPI fraction"
  ) +
  ylim(0.7,1) +
  theme_minimal()






####T cell exhaustion/suppression from extracellular lactate accumulation only:
###Sim func where we draw params from biologically informed/approximated distributions
simulate_profiles <- function(R, N = 400) {
  x <- seq(1e-4, R, length.out = N)
  dx <- diff(x)[1]
  
  #Basic geometry; how radius scales with surface area & vol.
  vol <- (4/3) * pi * R^3
  surf <- 4 * pi * R^2
  SA_V <- surf / vol  #surface area to volume ratio 
  
  #Params
  ##production/consumption per unit volume
  p_L0 <- 1     ##lactate produced per unit tumour volume
  p_T0 <- 1     ##O2 consumed per unit volume
  p_I0 <- 1     ##O2 consumed when immune cells are present at I = 1 (locally)
  
  ##boundary conditions scale with SA:Vol
  I_ext_base <- 1 ##functioning immune availability at tumour periphery, normalised to 1 to represent functioing immune presence at tumour edge
  O_ext_base <- 1 #normalised; assumed normoxic interface with surrounding env. (e.g, vasculature)
  L_ext_base <- pmax(0, 1 - SA_V) ##poorer lactate clearance with decreasing SA:Vol; captures the idea that larger tumours cannot clear lactate as efficiently, even at the boundary
  
  #diffusion const. (fixed)
  D_L <- 1 #normalised, baseline
  D_I <- 0.01 #immune cells infiltrate much slower/need active migration as they are bigger (constrained e.g., by stroma)
  D_O <- 10 #diffusion of O2 has to be at least one order of magnitude higher than those of larger species lactate or immune cells
  
  #Immune suppression by lactate ( Lactic acid inhibits T cell motility and cytokine release + induces transcriptional suppression of effector programs)
  #Lactate suppresses CTL activation; reversible by pH buffering
  a_L <- runif(1, 0.1, 0.9) #Captures immunosuppressive effects of low pH / lactic acidosis on T cells, via e.g., acid-sensitive signaling or T cell arrest (10-90% of normal)
  a_O <- 0#runif(1, 0.1, 0.9) #captures T cell exhaustion/suppression from (real) hypoxia (10-90% of normal activity)
  
  ##Thresholds for MPI
  O_thresh <- runif(1, 0.1, 0.5) #hypoxia is often defined within the range of 10%-50% of O2 conc. of normoxia (i.e., normal oxygen conc in healthy tissue)
  I_thresh <- 0.2 #Immune cell density below which tumour is considered "immune suppressed" for our score (MPI)
  ##so immune exhaustion when we have only 0.2 of the functioning immune cell to tumour cell saturation/equilibrium state (i..e, 1)
  
  ##Scaled params
  p_L <- p_L0 * vol ##lactate production scales with total metabolic demand, which grows with number of cells = volume
  p_T <- p_T0 * vol ##O2 consumption scales with total metabolic demand, which grows with number of cells = volume
  p_I <- p_I0 #fixed rate. While for tumour, total consumption scales with volume, for immune cells: localised consumption scales with local density
  I_ext <- I_ext_base * SA_V #Immune cells enter tumour from the outside (e.g., via lymph or vasculature), so immune access is also constrained by surface area to volume ratio
  O_ext <- O_ext_base * SA_V #Oxygen enters tumour via surface (e.g., blood vessels or capillary exchange). As volume increases faster than surface area, the same surface must serve a larger mass
  
  ##initialise profiles
  L <- rep(0, N) #assumes no pre-acidification
  I <- rep(0.1, N) #represents gradual infiltration or boundary supply
  O <- rep(1, N) #high initial value needed to allow depletion toward core
  
  #Looped solve
  tol <- 1e-4
  for (iter in 1:1000) {
    L_old <- L
    I_old <- I
    O_old <- O
    
    for (i in 2:(N-1)) {
      r <- x[i]
      
      #Lactate diffusion + production (discretised steady-state spherical diffusion eq.)
      L[i] <- ((r + dx/2)^2 * L[i+1] + (r - dx/2)^2 * L[i-1] + dx^2 * p_L * r^2 / D_L) /
        ((r + dx/2)^2 + (r - dx/2)^2)
      #^ lactate diffuses from production sources (tumour cells) 
      #^uniform production throughout the tumour (irregardless of local O2 conc), scaled by tumour volume
      #^lactate then diffuses outward but clearance is harder in larger tumours due to higher volume but falling SA
      
      ##Immune diffusion with suppression
      hypoxia_penalty <- ifelse(O[i] < O_thresh, 1, 0)
      denom <- ((r + dx/2)^2 + (r - dx/2)^2 + dx^2 * r^2 / D_I * (a_L * L[i] + a_O * hypoxia_penalty))
      #denom <- ((r + dx/2)^2 + (r - dx/2)^2 + dx^2 * a_L * L[i] * r^2 / D_I) #if not accounting for immune cell exhaustion
      #denom <- ((r + dx/2)^2 + (r - dx/2)^2 + dx^2 * r^2 / D_I * (a_L * L[i] + a_O * (1 - O[i])))  #if accounting for exhaustion smoothly
      I[i] <- ((r + dx/2)^2 * I[i+1] + (r - dx/2)^2 * I[i-1]) / denom
      #^immune cells diffuse inward from the tumour boundary (where we assume functioning vasculature)
      #^movement suppressed by high lactate conc. (acidosis), simulating T cell exhaustion/dysfunction
      #^if a_L is 0, this reduces to no lactate-mediated immune exhaustion (implicit suppression for more numerical stability)
      #^when O2 falls below threshold (O_thresh), hypoxia may occur, which can exhaust immune cells
      #Note that I[x] is a dimensionless/normalised concentration representing relative immune cell density at location x
      #It is not discrete, as we define it as a fractional occupancy (where 1 = maximum functioning immune infiltration [representing e.g., 1 immune cell per tumour cell], 0 = no immune infiltration)
      
      
      
      # Oxygen diffusion with consumption (oxygen consumption from immune cells, proportional to their local presence, and oxygen consumption cancer cells scaling with total volume)
      total_consump <- p_T + p_I * I[i] #oxygen demand from immune cells is local, and scaled by relative density
      O[i] <- ((r + dx/2)^2 * O[i+1] + (r - dx/2)^2 * O[i-1] - dx^2 * total_consump * r^2 / D_O) /
        ((r + dx/2)^2 + (r - dx/2)^2)
      ##O2 is supplied from the outside (tumour boundary/surface area), as we assume functional vasculature there
      #it diffuses inward but is consumed by both tumour cells and infiltrating immune cells
    }
    
    ##boundary conditions
    L[N] <- L_ext_base  ##if 0, we assume lactate is rapidly cleared at the tumour edge regardless of size
    L[1] <- L[2]
    I[N] <- I_ext
    I[1] <- I[2]
    O[N] <- O_ext
    O[1] <- O[2]
    
    if (max(abs(L - L_old), abs(I - I_old), abs(O - O_old)) < tol) break
  }
  
  ##MPI zone: when functioning immune cells cannot infiltrate accordingly/ infiltrated cells become exhausted/dysfunctional
  MPI_zone <- which(I < I_thresh)
  MPI_frac <- length(MPI_zone) / N
  
  tibble(
    x = x, L = L, I = I, O = O,
    MPI = x %in% x[MPI_zone],
    radius = R,
    MPI_fraction = MPI_frac
  )
}

#sim settings
radii <- seq(0.5, 5, by = 0.25)
n_reps <- 50  #param draws per tested radius

#run
set.seed(42)
all_results <- map_dfr(radii, function(r) {
  map_dfr(1:n_reps, function(rep) {
    simulate_profiles(R = r) %>%
      mutate(sim_id = rep)
  })
})

##conc profiles
results_long <- all_results %>%
  pivot_longer(cols = c(L, I, O), names_to = "Component", values_to = "Value")

ggplot(results_long, aes(x = x, y = Value, color = Component)) +
  geom_line(linewidth = 1) +
  facet_wrap(~radius, labeller = label_both, scales = "free_x") +
  labs(
    title = "Concentration profiles vs tumour radius",
    x = "Radial distance (cm)", y = "value"
  ) +
  theme_minimal()

##MPI emergence (not sure what else to call this)
mpi_summary <- all_results %>%
  group_by(radius) %>%
  summarise(
    median_MPI = median(MPI_fraction),
    lower = quantile(MPI_fraction, 0.25),
    upper = quantile(MPI_fraction, 0.75)
  )


ggplot(mpi_summary, aes(x = radius, y = median_MPI)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(
    name = "Tumour radius (cm)",
    sec.axis = sec_axis(~ . * 2, name = "Tumour diameter (cm)")
  ) +
  labs(
    title = "T cell exhaustion from lactate",
    y = "Median MPI fraction"
  ) +
  ylim(0.7,1) +
  theme_minimal()








####T cell exhaustion/suppression from extracellular lactate accumulation & true hypoxia combined:
###Sim func where we draw params from biologically informed/approximated distributions
simulate_profiles <- function(R, N = 400) {
  x <- seq(1e-4, R, length.out = N)
  dx <- diff(x)[1]
  
  #Basic geometry; how radius scales with surface area & vol.
  vol <- (4/3) * pi * R^3
  surf <- 4 * pi * R^2
  SA_V <- surf / vol  #surface area to volume ratio 
  
  #Params
  ##production/consumption per unit volume
  p_L0 <- 1     ##lactate produced per unit tumour volume
  p_T0 <- 1     ##O2 consumed per unit volume
  p_I0 <- 1     ##O2 consumed when immune cells are present at I = 1 (locally)
  
  ##boundary conditions scale with SA:Vol
  I_ext_base <- 1 ##functioning immune availability at tumour periphery, normalised to 1 to represent functioing immune presence at tumour edge
  O_ext_base <- 1 #normalised; assumed normoxic interface with surrounding env. (e.g, vasculature)
  L_ext_base <- pmax(0, 1 - SA_V) ##poorer lactate clearance with decreasing SA:Vol; captures the idea that larger tumours cannot clear lactate as efficiently, even at the boundary
  
  #diffusion const. (fixed)
  D_L <- 1 #normalised, baseline
  D_I <- 0.01 #immune cells infiltrate much slower/need active migration as they are bigger (constrained e.g., by stroma)
  D_O <- 10 #diffusion of O2 has to be at least one order of magnitude higher than those of larger species lactate or immune cells
  
  #Immune suppression by lactate ( Lactic acid inhibits T cell motility and cytokine release + induces transcriptional suppression of effector programs)
  #Lactate suppresses CTL activation; reversible by pH buffering
  a_L <- runif(1, 0.1, 0.9) #Captures immunosuppressive effects of low pH / lactic acidosis on T cells, via e.g., acid-sensitive signaling or T cell arrest (10-90% of normal)
  a_O <- runif(1, 0.1, 0.9) #captures T cell exhaustion/suppression from (real) hypoxia (10-90% of normal activity)
  
  ##Thresholds for MPI
  O_thresh <- runif(1, 0.1, 0.5) #hypoxia is often defined within the range of 10%-50% of O2 conc. of normoxia (i.e., normal oxygen conc in healthy tissue)
  I_thresh <- 0.2 #Immune cell density below which tumour is considered "immune suppressed" for our score (MPI)
  ##so immune exhaustion when we have only 0.2 of the functioning immune cell to tumour cell saturation/equilibrium state (i..e, 1)
  
  ##Scaled params
  p_L <- p_L0 * vol ##lactate production scales with total metabolic demand, which grows with number of cells = volume
  p_T <- p_T0 * vol ##O2 consumption scales with total metabolic demand, which grows with number of cells = volume
  p_I <- p_I0 #fixed rate. While for tumour, total consumption scales with volume, for immune cells: localised consumption scales with local density
  I_ext <- I_ext_base * SA_V #Immune cells enter tumour from the outside (e.g., via lymph or vasculature), so immune access is also constrained by surface area to volume ratio
  O_ext <- O_ext_base * SA_V #Oxygen enters tumour via surface (e.g., blood vessels or capillary exchange). As volume increases faster than surface area, the same surface must serve a larger mass
  
  ##initialise profiles
  L <- rep(0, N) #assumes no pre-acidification
  I <- rep(0.1, N) #represents gradual infiltration or boundary supply
  O <- rep(1, N) #high initial value needed to allow depletion toward core
  
  #Looped solve
  tol <- 1e-4
  for (iter in 1:1000) {
    L_old <- L
    I_old <- I
    O_old <- O
    
    for (i in 2:(N-1)) {
      r <- x[i]
      
      #Lactate diffusion + production (discretised steady-state spherical diffusion eq.)
      L[i] <- ((r + dx/2)^2 * L[i+1] + (r - dx/2)^2 * L[i-1] + dx^2 * p_L * r^2 / D_L) /
        ((r + dx/2)^2 + (r - dx/2)^2)
      #^ lactate diffuses from production sources (tumour cells) 
      #^uniform production throughout the tumour (irregardless of local O2 conc), scaled by tumour volume
      #^lactate then diffuses outward but clearance is harder in larger tumours due to higher volume but falling SA
      
      ##Immune diffusion with suppression
      hypoxia_penalty <- ifelse(O[i] < O_thresh, 1, 0)
      denom <- ((r + dx/2)^2 + (r - dx/2)^2 + dx^2 * r^2 / D_I * (a_L * L[i] + a_O * hypoxia_penalty))
      #denom <- ((r + dx/2)^2 + (r - dx/2)^2 + dx^2 * a_L * L[i] * r^2 / D_I) #if not accounting for immune cell exhaustion
      #denom <- ((r + dx/2)^2 + (r - dx/2)^2 + dx^2 * r^2 / D_I * (a_L * L[i] + a_O * (1 - O[i])))  #if accounting for exhaustion smoothly
      I[i] <- ((r + dx/2)^2 * I[i+1] + (r - dx/2)^2 * I[i-1]) / denom
      #^immune cells diffuse inward from the tumour boundary (where we assume functioning vasculature)
      #^movement suppressed by high lactate conc. (acidosis), simulating T cell exhaustion/dysfunction
      #^if a_L is 0, this reduces to no lactate-mediated immune exhaustion (implicit suppression for more numerical stability)
      #^when O2 falls below threshold (O_thresh), hypoxia may occur, which can exhaust immune cells
      #Note that I[x] is a dimensionless/normalised concentration representing relative immune cell density at location x
      #It is not discrete, as we define it as a fractional occupancy (where 1 = maximum functioning immune infiltration [representing e.g., 1 immune cell per tumour cell], 0 = no immune infiltration)
      
      
      
      # Oxygen diffusion with consumption (oxygen consumption from immune cells, proportional to their local presence, and oxygen consumption cancer cells scaling with total volume)
      total_consump <- p_T + p_I * I[i] #oxygen demand from immune cells is local, and scaled by relative density
      O[i] <- ((r + dx/2)^2 * O[i+1] + (r - dx/2)^2 * O[i-1] - dx^2 * total_consump * r^2 / D_O) /
        ((r + dx/2)^2 + (r - dx/2)^2)
      ##O2 is supplied from the outside (tumour boundary/surface area), as we assume functional vasculature there
      #it diffuses inward but is consumed by both tumour cells and infiltrating immune cells
    }
    
    ##boundary conditions
    L[N] <- L_ext_base  ##if 0, we assume lactate is rapidly cleared at the tumour edge regardless of size
    L[1] <- L[2]
    I[N] <- I_ext
    I[1] <- I[2]
    O[N] <- O_ext
    O[1] <- O[2]
    
    if (max(abs(L - L_old), abs(I - I_old), abs(O - O_old)) < tol) break
  }
  
  ##MPI zone: when functioning immune cells cannot infiltrate accordingly/ infiltrated cells become exhausted/dysfunctional
  MPI_zone <- which(I < I_thresh)
  MPI_frac <- length(MPI_zone) / N
  
  tibble(
    x = x, L = L, I = I, O = O,
    MPI = x %in% x[MPI_zone],
    radius = R,
    MPI_fraction = MPI_frac
  )
}

#sim settings
radii <- seq(0.5, 5, by = 0.25)
n_reps <- 50  #param draws per tested radius

#run
set.seed(42)
all_results <- map_dfr(radii, function(r) {
  map_dfr(1:n_reps, function(rep) {
    simulate_profiles(R = r) %>%
      mutate(sim_id = rep)
  })
})

##conc profiles
results_long <- all_results %>%
  pivot_longer(cols = c(L, I, O), names_to = "Component", values_to = "Value")

ggplot(results_long, aes(x = x, y = Value, color = Component)) +
  geom_line(linewidth = 1) +
  facet_wrap(~radius, labeller = label_both, scales = "free_x") +
  labs(
    title = "Concentration profiles vs tumour radius",
    x = "Radial distance (cm)", y = "value"
  ) +
  theme_minimal()

##MPI emergence (not sure what else to call this)
mpi_summary <- all_results %>%
  group_by(radius) %>%
  summarise(
    median_MPI = median(MPI_fraction),
    lower = quantile(MPI_fraction, 0.25),
    upper = quantile(MPI_fraction, 0.75)
  )


ggplot(mpi_summary, aes(x = radius, y = median_MPI)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(
    name = "Tumour radius (cm)",
    sec.axis = sec_axis(~ . * 2, name = "Tumour diameter (cm)")
  ) +
  labs(
    title = "T cell exhaustion from hypoxia & lactate",
    y = "Median MPI fraction"
  ) +
  ylim(0.7,1) +
  theme_minimal()












