# Matt Prill
# 14/04/24
# prill.mattfb@gmail.com
# Dissertation
# Data Analysis Code


# Libraries ----
library(tidyverse)  # Dplyr etc.
library(vegan)  # Diversity function
library(brms)  # Bayesian Modelling


getwd()
# Data ----
# Monthly Dagger Summary Stats
dagger_monthly <- read.csv("Data/daggers/imported_daggers/all/daggers_summary_monthly.csv")

dagger_summary_stats  <- read.csv("Data/daggers/imported_daggers/all/daggers_summary_stats.csv")


daggers_total <- read.csv("Data/daggers/imported_daggers/all/all_correct_times_4to1pm.csv") %>%
  mutate_at(vars(1:4), factor)


boghall <- read.csv("Data/invertebrates/boghall_inverts.csv",
                    skip = 1) %>%  # omit first row 
  mutate_at(vars(1:4), factor)  # Set appropriate columns as factors


adult_larvae_merged <- read.csv("Data/invertebrates/adult_larvae_merged.csv",
                                skip = 1) %>%  # omit first row 
  mutate_at(vars(1:4), factor)  # Set appropriate columns as factors


merged_abundances_only <- adult_larvae_merged %>% 
  select(-c(1:4)) %>%   # Removing non-data
  mutate_all(~coalesce(., 0)) # Replacing NAs with 0s

merged_abundances <- adult_larvae_merged %>% 
  mutate_all(~replace_na(., 0))
  

# Abundance data only
boghall_abundances <- boghall %>% 
  select(-c(1:4)) %>%   # Removing non-data
  mutate_all(~coalesce(., 0)) # Replacing NAs with 0s

# relative abundances (carabids)
relative_abundances <- read.csv("Data/invertebrates/relative_abundances.csv")
relative_abundances <- relative_abundances %>% 
  mutate_all(~coalesce(., 0))


# Surface Pitfalls
above_ground <- boghall %>% 
  filter(grepl("^pitfall_", trap))


# Subterranean Traps
below_ground <- boghall %>% 
  filter(grepl("^sub_", trap))


# Above ground (Carabids)
above_ground_carabids <- above_ground %>% 
  mutate(trap_total = rowSums(select(., 18:43), na.rm = TRUE)) %>%  # Summing each row carabids (adult and larvae)
  select(1:4, 18:43, 49)

# abundances
above_ground_abundances <- above_ground_carabids %>% 
  select(5:30) %>% 
  mutate_all(~replace(., is.na(.), 0))  # replace na with 0s




# Below ground (Carabids)
below_ground_carabids <- below_ground %>% 
  mutate(trap_total = rowSums(select(., 18:43), na.rm = TRUE)) %>%  # Summing each row carabids (adult and larvae)
  select(1:4, 18:43, 49)

below_ground_abundances <- below_ground_carabids %>% 
  select(5:30) %>% 
  mutate_all(~replace(., is.na(.), 0))  # replace na with 0



# ALL SLUGS
slugs <- boghall %>% 
  select(1:4, 8) %>% 
  mutate_all(~replace(., is.na(.), 0))  # replace na with 0s


# Above ground Slugs
above_ground_slugs <- above_ground %>% 
  select(1:4, 8)

# Below ground slugs
below_ground_slugs <- below_ground %>% 
  select(1:4,8)


# Raw data
adult_larvae_merged <- read.csv("Data/invertebrates/adult_larvae_merged.csv",
                                skip = 1) %>%  # omit first row 
  mutate_at(vars(1:4), factor)  # Set appropriate columns as factors

# Above and Below
above_merged <- adult_larvae_merged %>% 
  filter(grepl("^pitfall_", trap))

below_merged <-adult_larvae_merged %>% 
  filter(grepl("^sub_", trap))



above_merged_abundances_only <- above_merged %>% 
  select(-c(1:4)) %>%   # Removing non-data
  mutate_all(~coalesce(., 0)) # Replacing NAs with 0s

below_merged_abundances_only <- below_merged %>% 
  select(-c(1:4)) %>%   # Removing non-data
  mutate_all(~coalesce(., 0)) # Replacing NAs with 0s

above_merged_abundances_only_months <- above_merged %>% 
  mutate_all(~replace_na(., 0))

below_merged_abundances_only_months <- below_merged %>% 
  mutate_all(~replace_na(., 0))

# Relative abundances
above_merged_relative_abundances <- decostand(above_merged_abundances_only, method = "hellinger")

below_merged_relative_abundances <- decostand(above_merged_abundances_only, method = "hellinger")


# Richness, diverstity and total abundance 

above_merged_divs <- above_merged_abundances_only_months %>%
  mutate(shannons_diversity = diversity(above_merged_abundances_only),  # Add column of shannons div
         richness = specnumber(above_merged_abundances_only),  # Add column of species richness
         total_abundance = rowSums(above_merged_abundances_only)) %>%   # Total abundance
  select(month, plot_type, plot, trap, shannons_diversity, richness, total_abundance)  # Remove abundance data


below_merged_divs <- below_merged_abundances_only_months %>%
  mutate(shannons_diversity = diversity(below_merged_abundances_only),  # Add column of shannons div
         richness = specnumber(below_merged_abundances_only),  # Add column of species richness
         total_abundance = rowSums(below_merged_abundances_only)) %>%   # Total abundance
  select(month, plot_type, plot, trap, shannons_diversity, richness, total_abundance)  # Remove abundance data





# ---- 
# ----
# ----
# Above ground carabid abundance Bayesian Analyses ----

# Data distribution - Histogram & Distribution Curve
(above_abundance_dist <- ggplot(above_ground_carabids, aes(x = trap_total)) +
   geom_histogram(binwidth = 5, color = "darkblue", fill="lightblue") +
   ylab("Count\n") +
   xlab("\n\n Pitfall carabid total ") +
   theme_classic() + # theme_bw() to add 
   theme(axis.text=element_text(size = 20),
         axis.title=element_text(size = 20))) +
  geom_density(aes(y=8.1*..count..))


# Data Distribution, (histogram and distribution curve (latter grouped by plot type))
(above_abundance_dist_both <- ggplot(above_ground_carabids, aes(x = trap_total, fill = plot_type)) +
    geom_histogram(binwidth = 5, color = "darkblue", position = "identity", alpha = 0.5) +
    ylab("Count\n") +
    xlab("\n\n Pitfall carabid total ") +
    theme_classic() + # theme_bw() to add 
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20)) +
    geom_density(aes(y = 8.1 * ..count.., color = plot_type), alpha = 0.5) + # distribution curve
    scale_fill_manual(values = c("lightblue", "lightgreen")) + 
    scale_color_manual(values = c("blue", "green")))

# Normality Test
shapiro.test(above_ground_carabids$trap_total)  # Doesn't need to be normal

# Bayesian analysis 
above_abundance_mbrms <- brm(trap_total ~ plot_type + month + (1|plot),  # pop = response variable, plot_type + month = explanatory variables, plot = random effect
                             data = above_ground_carabids, family = poisson(), chains = 3,  # Family argument must reflect the data distribution
                             iter = 3000, warmup = 1000)


# Assessing model fit
plot(above_abundance_mbrms)  # Plot the model
pp_check(above_abundance_mbrms) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(above_abundance_mbrms)  # Summaries model outputs


# The summary will be different each time the model is run due to stochasticity 
# Under 'population-level effects' the 'Estimate' gives the mean of the posterior distribution for each variable 
# They can be used as the intercept and slope for the relationship between our two variables.
# Est.Error is the error associated with those means (the standard error).
# Credibility Interval (CI) = The interval where 95% of posterior distribution values fall
# KEY! if the Interval encompasses 0, the the effect is not significant, if it is strictly positive or negative, it is significant 
# Bulk_EES and Tail_EES are the sample size measures for each parameter and should be > 1000
# Rhat should be 1 for each effect if the model has converged well


# Back transforming results
# We need to here because the model log-transformed the data because of the poisson distribution.
# Add the mean to the intercept output = actual value of increase
# 0.32 + 2.3 = 2.62
# Get the exponential of that value to undo the log-transformation
# exp(2.62)
# = 13.73572
exp(2.62)

m <- subset(above_ground_carabids, plot_type == "mustard")
f <- subset(above_ground_carabids, plot_type == "fallow")
mean(m$trap_total)
mean(f$trap_total)
max(m$trap_total)
max(f$trap_total)

# bayesian analysis (NO RANDOM EFFECT)
above_abundance_mbrms_2 <- brm(trap_total ~ plot_type + month,  # pop = response variable, plot_type + month = explanatory variables, plot = random effect
                               data = above_ground_carabids, family = poisson(), chains = 3,  # Family argument must reflect the data distribution
                               iter = 3000, warmup = 1000)


# Assessing model fit
plot(above_abundance_mbrms_2)  # Plot the model
pp_check(above_abundance_mbrms_2) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(above_abundance_mbrms_2)  # Summaries model outputs



# bayesian analysis 
above_abundance_mbrms_3 <- brm(trap_total ~ plot_type + month + (1|plot/trap),  # pop = response variable, plot_type + month = explanatory variables, plot = random effect
                               data = above_ground_carabids, family = poisson(), chains = 3,  # Family argument must reflect the data distribution
                               iter = 3000, warmup = 1000)


# Assessing model fit
plot(above_abundance_mbrms_3)  # Plot the model
pp_check(above_abundance_mbrms_3) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(above_abundance_mbrms_3)  # Summaries model outputs


# LOO validation - ASSESSING + COMPARING THE FIT OF THE MODELS 
# Based on elpd (higher value = better fit)
loo(above_abundance_mbrms, above_abundance_mbrms_2, above_abundance_mbrms_3, compare = TRUE)  # First model is best fit (higher elpd_diff)


# Below ground carabid abundance Bayesian Analyses ----

# Data Distribution Visualisation (Clear Poission)
(below_abundance_dist <- ggplot(below_ground_carabids, aes(x = trap_total)) +
   geom_histogram(binwidth = 1, color = "darkblue", fill="lightblue") +
   ylab("Count\n") +
   xlab("\n Pitfall carabid total ") +
   theme_classic() + # theme_bw() to add 
   theme(axis.text=element_text(size = 20),
         axis.title=element_text(size = 20)))


# If I suspect that the relationship between plot type and abundance might differ,
# it could be beneficial to visualize the data separately for each plot type (not the case):

# Data of distribution by plot type
(below_abundance_dist_both <- ggplot(below_ground_carabids, aes(x = trap_total, fill = plot_type)) +
    geom_histogram( color = "darkblue", position = "identity", alpha = 0.5) +
    ylab("Count\n") +
    xlab("\n\n Pitfall carabid total ") +
    theme_classic() + # theme_bw() to add 
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20)) +
    geom_density(aes(y = 1 * ..count.., color = plot_type), alpha = 0.5) + # distribution curve
    scale_fill_manual(values = c("lightblue", "lightgreen")) + # Choose appropriate colors
    scale_color_manual(values = c("blue", "green"))) # Choose appropriate colors for density plot


# Normality Test
shapiro.test(below_ground_carabids$trap_total)  # Doesn't need to be normal


# Bayesian analysis
below_abundance_mbrms <- brm(trap_total ~ plot_type + month + (1|plot),  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                             data = below_ground_carabids, family = poisson(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                             iter = 3000, warmup = 1000)

# Assessing model fit
plot(below_abundance_mbrms)  # Plot the model
pp_check(below_abundance_mbrms) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(below_abundance_mbrms)  # Summaries model outputs



# Back transforming results
# We need to here because the model log-transformed the data because of the poisson distribution.
# Add the mean to the intercept output = actual value of increase
0.28 + 1.27
# Get the exponential of that value to undo the log-transformation
exp(1.55)
# = 4.706648




# Bayesian analysis (NO RANDOM EFFECT)
below_abundance_mbrms_2 <- brm(trap_total ~ plot_type + month,  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                               data = below_ground_carabids, family = poisson(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                               iter = 3000, warmup = 1000)

# Assessing model fit
plot(below_abundance_mbrms_2)  # Plot the model
pp_check(below_abundance_mbrms_2) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(below_abundance_mbrms_2)  # Summaries model outputs



# IGNORE THIS
# Back transforming results
# We need to here because the model log-transformed the data because of the poisson distribution.
# Add the mean to the intercept output = actual value of increase
# 0.29 + 1.28 = 1.57
# Get the exponential of that value to undo the log-transformation
# exp(1.57)
# = 4.806648
exp(1.57) # IGNORE


below_abundance_mbrms_3 <- brm(trap_total ~ plot_type + month + (1|plot/trap),  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                               data = below_ground_carabids, family = poisson(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                               iter = 3000, warmup = 1000)

# Assessing model fit
plot(below_abundance_mbrms_3)  # Plot the model
pp_check(below_abundance_mbrms_3) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(below_abundance_mbrms_3)  # Summaries model outputs



# LOO validation - ASSESSING + COMPARING THE FIT OF THE MODELS 
# Based on elpd (higher value = better fit)
loo(below_abundance_mbrms, below_abundance_mbrms_2, below_abundance_mbrms_3, compare = TRUE)  # 3 is better fit but insufficient no of 'trap' levels





# ----
# ----
# ----
# Pest Analyses ----

# Visualising distribution, poisson
(above_abundance_dist <- ggplot(slugs, aes(x = stylommatophora)) +
   geom_histogram(color = "darkblue", fill="lightblue") +
   ylab("Count\n") +
   xlab("\n\n Pitfall carabid total ") +
   theme_classic() + # theme_bw() to add 
   theme(axis.text=element_text(size = 20),
         axis.title=element_text(size = 20))) +
  geom_density(aes(y=1*..count..)) # distribution curve

# bayesian analysis 
slug_abundance_mbrms <- brm(stylommatophora ~ plot_type + month,  # pop = response variable, plot_type + month = explanatory variables, plot = random effect
                            data = slugs, family = poisson(), chains = 3,  # Family argument must reflect the data distribution
                            iter = 3000, warmup = 1000)
# back transformation
0.8 + 0.03 
exp(0.83)
# 2.29 more slugs on avg in fallow plots


# Assessing model fit
plot(slug_abundance_mbrms)  # Plot the model
pp_check(slug_abundance_mbrms) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(slug_abundance_mbrms)  # Summaries model outputs


slug_abundance_mbrms_2 <- brm(stylommatophora ~ plot_type + month + (1|plot),  # pop = response variable, plot_type + month = explanatory variables, plot = random effect
                              data = slugs, family = poisson(), chains = 3,  # Family argument must reflect the data distribution
                              iter = 3000, warmup = 1000)

# Assessing model fit
plot(slug_abundance_mbrms_2)  # Plot the model
pp_check(slug_abundance_mbrms_2) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(slug_abundance_mbrms_2)  # Summaries model outputs


slug_abundance_mbrms_3 <- brm(stylommatophora ~ plot_type + month + (1|plot/trap),  # pop = response variable, plot_type + month = explanatory variables, plot = random effect
                              data = slugs, family = poisson(), chains = 3,  # Family argument must reflect the data distribution
                              iter = 3000, warmup = 1000)

# Assessing model fit
plot(slug_abundance_mbrms_3)  # Plot the model
pp_check(slug_abundance_mbrms_3) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(slug_abundance_mbrms_3)  # Summaries model outputs


# LOO validation - ASSESSING + COMPARING THE FIT OF THE MODELS 
# Based on elpd (higher value = better fit)
loo(slug_abundance_mbrms, slug_abundance_mbrms_2, slug_abundance_mbrms_3, compare = TRUE)  # Second model is best fit (higher elpd_diff and pp fit better, no divergent transitions)
loo(slug_abundance_mbrms, slug_abundance_mbrms_2, slug_abundance_mbrms_3)  # Second model is best fit (higher elpd_diff and pp fit better, no divergent transitions)


# ----
# ----
# ----
# Above Richness Analysis ----
# Distribution
(above_richness_dist <- ggplot(above_merged_divs, aes(x= richness)) +
   geom_histogram(binwidth = 1, color = "darkblue", fill="lightblue") +
   ylab("Count\n") +
   xlab("\n\n Pitfall carabid total ") +
   theme_classic() + # theme_bw() to add 
   theme(axis.text=element_text(size = 20),
         axis.title=element_text(size = 20))) +
  geom_density(aes(y=1*..count..)) # distribution curve



# Bayesian Anlaysis
above__merged_richness_mbrms <- brm(richness ~ plot_type + month + (1|plot),  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                                    data = above_merged_divs, family = gaussian(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                                    iter = 3000, warmup = 1000)

plot(above__merged_richness_mbrms)  # Plot the model
pp_check(above__merged_richness_mbrms) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(above__merged_richness_mbrms)  # Summaries model outputs


# Bayesian Anlaysis (NO RANDOM EFFECT)
above_merged_richness_mbrms_2 <- brm(richness ~ plot_type + month,  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                                     data = above_merged_divs, family = gaussian(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                                     iter = 3000, warmup = 1000)

plot(above_merged_richness_mbrms_2)  # Plot the model
pp_check(above_merged_richness_mbrms_2) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(above_merged_richness_mbrms_2)  # Summaries model outputs


# Bayesian Anlaysis (NO RANDOM EFFECT)
above_merged_richness_mbrms_3 <- brm(richness ~ plot_type + month + (1|plot/trap),  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                                     data = above_merged_divs, family = gaussian(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                                     iter = 3000, warmup = 1000)

plot(above_merged_richness_mbrms_3)  # Plot the model
pp_check(above_merged_richness_mbrms_3) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(above_merged_richness_mbrms_3)  # Summaries model outputs

# LOO validation - ASSESSING + COMPARING THE FIT OF THE MODELS 
# Based on elpd (higher value = better fit)
loo(above_merged_richness_mbrms, above_merged_richness_mbrms_2, above_merged_richness_mbrms_3, compare = TRUE)  # Second model is best fit (higher elpd_diff and pp fit better, no divergent transitions)


sum(above_merged_divs$plot_type == "mustard")
sum(above_merged_divs$plot_type == "fallow")

x = subset(above_merged_divs, plot_type == "mustard")
sd(x$richness)
sd(x$shannons_diversity)

y = subset(above_merged_divs, plot_type == "fallow")
sd(y$richness)
sd(y$shannons_diversity)


x = subset(below_merged_divs, plot_type == "mustard")
sd(x$richness)
sd(x$shannons_diversity)

y = subset(below_merged_divs, plot_type == "fallow")
sd(y$richness)
sd(y$shannons_diversity)

# Below Richness Analysis ----
# Distribution
(above_richness_dist <- ggplot(below_merged_divs, aes(x= richness)) +
   geom_histogram(binwidth = 1, color = "darkblue", fill="lightblue") +
   ylab("Count\n") +
   xlab("\n\n Pitfall carabid total ") +
   theme_classic() + # theme_bw() to add 
   theme(axis.text=element_text(size = 20),
         axis.title=element_text(size = 20))) +
  geom_density(aes(y=1*..count..)) # distribution curve
below__merged_richness_mbrms <- brm(richness ~ plot_type + month + (1|plot),  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                                    data = below_merged_divs, family = gaussian(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                                    iter = 3000, warmup = 1000)

plot(below__merged_richness_mbrms)  # Plot the model
pp_check(below__merged_richness_mbrms) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(below__merged_richness_mbrms)  # Summaries model outputs


# Bayesian Anlaysis (NO RANDOM EFFECT)
below__merged_richness_mbrms_2 <- brm(richness ~ plot_type + month,  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                                      data = below_merged_divs, family = gaussian(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                                      iter = 3000, warmup = 1000)

plot(below__merged_richness_mbrms_2)  # Plot the model
pp_check(below__merged_richness_mbrms_2) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(below__merged_richness_mbrms_2)  # Summaries model outputs


# Bayesian Anlaysis (NO RANDOM EFFECT)
below__merged_richness_mbrms_3 <- brm(richness ~ plot_type + month + (1|plot/trap),  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                                      data = below_merged_divs, family = gaussian(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                                      iter = 3000, warmup = 1000)

plot(below__merged_richness_mbrms_3)  # Plot the model
pp_check(below__merged_richness_mbrms_3) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(below__merged_richness_mbrms_3)  # Summaries model outputs

# LOO validation - ASSESSING + COMPARING THE FIT OF THE MODELS 
# Based on elpd (higher value = better fit)
loo(below__merged_richness_mbrms, below__merged_richness_mbrms_2, below__merged_richness_mbrms_3, compare = TRUE)  # Second model is best fit (higher elpd_diff and pp fit better, no divergent transitions)

# ----
# ----
# ----
# Above Shannons Diversity Analysis ----
# Ensure Gaussian distribution
(above_shannons_dist <- ggplot(above_merged_divs, aes(x = shannons_diversity)) +
   geom_histogram(color = "darkblue", fill="lightblue") +
   ylab("Count\n") +
   xlab("\n\n Shannons Div ") +
   #scale_y_continuous(expand = c(0,0), limits = c(0,7), breaks = seq(0,7,1)) +
   #scale_x_continuous(expand = c(0,0), limits = c(0,200), breaks = seq(0,200,25)) +
   theme_classic() + # theme_bw() to add 
   theme(axis.text=element_text(size = 20),
         axis.title=element_text(size = 20))) +
  geom_density(aes(y=.1*..count..)) # distribution curve


# Bayesian Analysis
above_merged_divs_mbrms <- brm(shannons_diversity ~ plot_type + month + (1|plot),  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                               data = above_merged_divs, family = gaussian(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                               iter = 3000, warmup = 1000)

plot(above_merged_divs_mbrms)  # Plot the model
pp_check(above_merged_divs_mbrms) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(above_merged_divs_mbrms)  # Summaries model outputs


# Bayesian Analysis (NO RANDOM EFFECT)
above_merged_divs_mbrms_2 <- brm(shannons_diversity ~ plot_type + month,  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                                 data = above_merged_divs, family = gaussian(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                                 iter = 3000, warmup = 1000)

plot(above_merged_divs_mbrms_2)  # Plot the model
pp_check(above_merged_divs_mbrms_2) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(above_merged_divs_mbrms_2)  # Summaries model outputs


# Bayesian Analysis nested
above_merged_divs_mbrms_3 <- brm(shannons_diversity ~ plot_type + month + (1|plot/trap),  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                                 data = above_merged_divs, family = gaussian(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                                 iter = 3000, warmup = 1000)

plot(above_merged_divs_mbrms_3)  # Plot the model
pp_check(above_merged_divs_mbrms_3) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(above_merged_divs_mbrms_3)  # Summaries model outputs


# LOO validation - ASSESSING + COMPARING THE FIT OF THE MODELS 
# Based on elpd (higher value = better fit)
loo(above_merged_divs_mbrms, above_merged_divs_mbrms_2, above_merged_divs_mbrms_3, compare = TRUE)  # Second model is best fit (higher elpd_diff and pp fit better, no divergent transitions)

# Below Shannons Diversity Analysis ----
(above_shannons_dist <- ggplot(below_merged_divs, aes(x = shannons_diversity)) +
   geom_histogram( color = "darkblue", fill="lightblue") +
   ylab("Count\n") +
   xlab("\n\n Shannons Div ") +
   #scale_y_continuous(expand = c(0,0), limits = c(0,7), breaks = seq(0,7,1)) +
   #scale_x_continuous(expand = c(0,0), limits = c(0,200), breaks = seq(0,200,25)) +
   theme_classic() + # theme_bw() to add 
   theme(axis.text=element_text(size = 20),
         axis.title=element_text(size = 20))) +
  geom_density(aes(y=.1*..count..)) # distribution curve


# Bayesian Analysis
below_merged_divs_mbrms <- brm(shannons_diversity ~ plot_type + month + (1|plot),  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                               data = below_merged_divs, family = gaussian(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                               iter = 3000, warmup = 1000)

plot(below_merged_divs_mbrms)  # Plot the model
pp_check(below_merged_divs_mbrms) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(below_merged_divs_mbrms)  # Summaries model outputs


# Bayesian Analysis (NO RANDOM EFFECT)
below_merged_divs_mbrms_2 <- brm(shannons_diversity ~ plot_type + month,  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                                 data = below_merged_divs, family = gaussian(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                                 iter = 3000, warmup = 1000)

plot(below_merged_divs_mbrms_2)  # Plot the model
pp_check(below_merged_divs_mbrms_2) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(below_merged_divs_mbrms_2)  # Summaries model outputs


# Bayesian Analysis nested
below_merged_divs_mbrms_3 <- brm(shannons_diversity ~ plot_type + month + (1|plot/trap),  # pop = response variable, plot_type = explanatory variable, plot = random effect. Month cant be random as it is explanatory and only has 3 levels (Sep,Oct,Nov)
                                 data = below_merged_divs, family = gaussian(), chains = 3,  # GONE FOR GAUSSIAN BC POISSION REQUIRES INTEGERS AND IS USED FOR COUNT DATA. LOOK INTO THE DISTRIBUTION OF THE SHANNONS DATA TO MAKE SURE MODEL IS GOOD
                                 iter = 3000, warmup = 1000)

plot(below_merged_divs_mbrms_3)  # Plot the model
pp_check(below_merged_divs_mbrms_3) # Showing posterior distribution (y) and 10 random 10 random distributions created by the model (yrep) 
summary(below_merged_divs_mbrms_3)  # Summaries model outputs


# LOO validation - ASSESSING + COMPARING THE FIT OF THE MODELS 
# Based on elpd (higher value = better fit)
loo(below_merged_divs_mbrms, below_merged_divs_mbrms_2, below_merged_divs_mbrms_3, compare = TRUE)  # Second model is best fit (higher elpd_diff and pp fit better, no divergent transitions)

# ----
# ----
# ----

# Above NMDS ----
above_merged_abundances_only

which(rowSums(above_merged_abundances_only) == 0)  # Row 53 = 0 (november, P2, mustard, pitfall_2)

above_merged_abundances_only_nmds <- above_merged_abundances_only %>% 
  filter(rowSums(.) > 0)  # Remove row that sums to 0 (ROW 53, 	november, P2 ,mustard, pitfall_2)

above_merged_abundances_only_nmds_rows <- above_merged_abundances_only_months %>% 
  slice(-c(53)) # Remove row 53 to make it mathc rda


# NMDS + fit checks
set.seed(123)
above_merged_nmds <- metaMDS(above_merged_abundances_only_nmds, distance = "bray", k = 3, autotransform = TRUE, trymax = 300, group = above_merged_abundances_only_months[, c("Plot_Type", "Month")]) # Bray as classic Euclidean distances are sensitive to species/order abundances, k = 3 for better fit
above_merged_nmds
plot(above_merged_nmds)
stressplot(above_merged_nmds)
above_merged_nmds$stress  # Stress is < 0.2 which is good

# Checking assumption of homogeneity in multivariate dispersion
above.inv.dist <- vegdist(above_merged_abundances_only_nmds, method = "bray")  
above.inv.dispersion <- betadisper(above.inv.dist, group = above_merged_abundances_only_nmds_rows$plot_type_month)
permutest(above.inv.dispersion)



# Analysis
above_merged_nmds_plot_type <- above_merged_abundances_only_nmds_rows$plot_type  # For analysis
above_merged_nmds_month <- above_merged_abundances_only_nmds_rows$month  # For analysis

above_merged_permanova <- adonis2(as.matrix(above_merged_abundances_only_nmds) ~ above_merged_nmds_plot_type + above_merged_nmds_month,
                                  permutations = 999, method = "bray")

above_merged_permanova  # Summary 




# Plotting

data.scores.above.merged <- as.data.frame(scores(above_merged_nmds)$sites)  # Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores.above.merged$plot_type <- as.factor(above_merged_abundances_only_nmds_rows$plot_type)  # Trying to group better
data.scores.above.merged$month <- as.factor(above_merged_abundances_only_nmds_rows$month)  # Trying to group better
data.scores.above.merged$site <- rownames(data.scores.above.merged)  # create a column of site names, from the rownames of data.scores
data.scores.above.merged$grp <- as.factor(above_merged_abundances_only_nmds_rows$plot_type)  #  add the group variable created earlier

data.scores.above.merged <- data.scores.above.merged %>%   # Multiple grouping variables
  mutate(plot_type_month = paste0(plot_type, "-", month))

head(data.scores.above.merged)

above_merged_abundances_only_nmds_rows <- above_merged_abundances_only_nmds_rows %>%         # Same for this data
  mutate(plot_type_month = paste0(plot_type, "-", month))


hull.data.above.merged <- data.frame()
for(i in 1:length(unique(above_merged_abundances_only_nmds_rows$plot_type_month))){
  
  temp.above.merged <- data.scores.above.merged[data.scores.above.merged$plot_type_month == unique(above_merged_abundances_only_nmds_rows$plot_type_month)[i], ][chull(data.scores.above.merged[data.scores.above.merged$plot_type_month == 
                                                                                                                                                                                                  unique(above_merged_abundances_only_nmds_rows$plot_type_month)[i], c("NMDS1", "NMDS2")]), ]
  hull.data.above.merged <- rbind(hull.data.above.merged, temp.above.merged)
}

# PLot

(invert_NMDS_plot <- ggplot() +
    geom_polygon(data=hull.data.above.merged, aes(x =NMDS1,y=NMDS2, fill = plot_type_month, group = plot_type_month), alpha = 0.50) + # add the convex hulls
    geom_point(data=data.scores.above.merged,aes(x = NMDS1, y = NMDS2, colour = plot_type_month), size = 8)) + # add the point markers
  scale_colour_manual("Plot Type - Month Combination", values = c("#63B0EBFE",  "slateblue2","darkblue","#6E4318", "saddlebrown", "#520E19FE")) +
  scale_fill_manual("Plot Type - Month Combination", values = c("#63B0EBFE", "slateblue2", "darkblue","#6E4318", "saddlebrown", "#520E19FE")) +
  theme_classic() +
  theme(axis.text = element_text(size = 50, colour = "black"),
        axis.text.x = element_text(size = 50, color = "black"),
        axis.title = element_text(size = 60),
        axis.ticks = element_line(size = 4, "black"),  # Adjusting axis ticks size
        axis.ticks.length = unit(0.25, "cm"),
        axis.line = element_line(size = 2),   # Adjusting axis line thickness
        panel.grid.major = element_blank(),  # Removing major grid lines
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 35),
        legend.title = element_text(size = 35)) +
  guides(fill = guide_legend(title = "Plot Type Month Combination",   # Change the legend title
                             keywidth = unit(3, "cm"),  # Set the width of the legend keys
                             keyheight = unit(2, "cm"),  # Set the height of the legend keys
                             title.position = "top",  # Position the title at the top
                             title.hjust = 0.5,  # Center-align the title
                             label.theme = element_text(size = 40)),  # Increase legend label size
         color = "none")  # Hide the color legend




# Below NMDS ----
below_merged_abundances_only

which(rowSums(below_merged_abundances_only) == 0)  # Rows 6, 9, 15, 38 = 0 (Sep, mustard x 3, nov, mustard )

below_merged_abundances_only_nmds <- below_merged_abundances_only %>% 
  filter(rowSums(.) > 0)  # Remove row that sums to 0 (ROW 53, 	november, P2 ,mustard, pitfall_2)

below_merged_abundances_only_nmds_rows <- below_merged_abundances_only_months %>% 
  slice(-c(6, 9, 15, 38 )) # Remove rows to make it match rda

length(below_merged_abundances_only_nmds_rows$plot_type)


# NMDS + fit checks
set.seed(123)
below_merged_nmds <- metaMDS(below_merged_abundances_only_nmds, distance = "bray", k = 3, autotransform = TRUE, trymax = 300, group = below_merged_abundances_only_months[, c("Plot_Type", "Month")]) # Bray as classic Euclidean distances are sensitive to species/order abundances, k = 3 for better fit
below_merged_nmds
plot(below_merged_nmds)
stressplot(below_merged_nmds)
below_merged_nmds$stress  # Stress is < 0.2 which is good


# Checking assumption of homogeneity in multivariate dispersion
below.inv.dist <- vegdist(below_merged_abundances_only_nmds, method = "bray") # 
below.inv.dispersion <- betadisper(below.inv.dist, group = below_merged_abundances_only_nmds_rows$plot_type_month)
permutest(below.inv.dispersion)



# Analysis
below_merged_nmds_plot_type <- below_merged_abundances_only_nmds_rows$plot_type  # For analysis
below_merged_nmds_month <- below_merged_abundances_only_nmds_rows$month  # For analysis

below_merged_permanova <- adonis2(as.matrix(below_merged_abundances_only_nmds) ~ below_merged_nmds_plot_type + below_merged_nmds_month, below_merged_abundances_only_nmds,
                                  permutations = 999, method = "bray")

below_merged_permanova  # Summary 




# Plotting

data.scores.below.merged <- as.data.frame(scores(below_merged_nmds)$sites)  # Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores.below.merged$plot_type <- as.factor(below_merged_abundances_only_nmds_rows$plot_type)  # Trying to group better
data.scores.below.merged$month <- as.factor(below_merged_abundances_only_nmds_rows$month)  # Trying to group better
data.scores.below.merged$site <- rownames(data.scores.below.merged)  # create a column of site names, from the rownames of data.scores
data.scores.below.merged$grp <- as.factor(below_merged_abundances_only_nmds_rows$plot_type)  #  add the group variable created earlier

data.scores.below.merged <- data.scores.below.merged %>%   # Multiple grouping variables
  mutate(plot_type_month = paste0(plot_type, "-", month))

head(data.scores.below.merged)

below_merged_abundances_only_nmds_rows <- below_merged_abundances_only_nmds_rows %>%         # Same for this data
  mutate(plot_type_month = paste0(plot_type, "-", month))


hull.data.below.merged <- data.frame()
for(i in 1:length(unique(below_merged_abundances_only_nmds_rows$plot_type_month))){
  
  temp.below.merged <- data.scores.below.merged[data.scores.below.merged$plot_type_month == unique(below_merged_abundances_only_nmds_rows$plot_type_month)[i], ][chull(data.scores.below.merged[data.scores.below.merged$plot_type_month == 
                                                                                                                                                                                                  unique(below_merged_abundances_only_nmds_rows$plot_type_month)[i], c("NMDS1", "NMDS2")]), ]
  hull.data.below.merged <- rbind(hull.data.below.merged, temp.below.merged)
}


(invert_NMDS_plot <- ggplot() +
    geom_polygon(data=hull.data.below.merged, aes(x =NMDS1,y=NMDS2, fill = plot_type_month, group = plot_type_month), alpha = 0.50) + # add the convex hulls
    geom_point(data=data.scores.below.merged,aes(x = NMDS1, y = NMDS2, colour = plot_type_month), size = 8)) + # add the point markers
  scale_colour_manual("Plot Type - Month Combination", values = c("#63B0EBFE",  "slateblue2","darkblue","#6E4318", "saddlebrown", "#520E19FE")) +
  scale_fill_manual("Plot Type - Month Combination", values = c("#63B0EBFE", "slateblue2", "darkblue","#6E4318", "saddlebrown", "#520E19FE")) +
  theme_classic() +
  theme(axis.text = element_text(size = 50, colour = "black"),
        axis.text.x = element_text(size = 50, color = "black"),
        axis.title = element_text(size = 60),
        axis.ticks = element_line(size = 4, "black"),  # Adjusting axis ticks size
        axis.ticks.length = unit(0.25, "cm"),
        axis.line = element_line(size = 2),   # Adjusting axis line thickness
        panel.grid.major = element_blank(),  # Removing major grid lines
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 35),
        legend.title = element_text(size = 35)) +
  guides(fill = guide_legend(title = "Plot Type Month Combination",   # Change the legend title
                             keywidth = unit(3, "cm"),  # Set the width of the legend keys
                             keyheight = unit(2, "cm"),  # Set the height of the legend keys
                             title.position = "top",  # Position the title at the top
                             title.hjust = 0.5,  # Center-align the title
                             label.theme = element_text(size = 40)),  # Increase legend label size
         color = "none")  # Hide the color legend

# ----
# ----
# ----
# Above Community Evenness ----
# Mustard Subsetting
above_merged_m <- subset(above_merged, plot_type == "mustard")

above_merged_relabun_m <- above_merged_m %>%  
  select(where(function(above_merged_m) any(!is.na(above_merged_m) & is.numeric(above_merged_m)))) %>% 
  mutate_all(~replace_na(., 0)) %>%             # Remove Na
  decostand(method = "hellinger") %>%           # Relative abund
  colMeans() %>%                                # Mean relative abundances
  as.data.frame() %>% 
  mutate(rank = rank(-.)) %>% 
  setNames(c("ra", "rank"))



# Fallow Subsetting
above_merged_f <- subset(above_merged, plot_type == "fallow")

above_merged_relabun_f <- above_merged_f %>%  
  select(where(function(above_merged_m) any(!is.na(above_merged_m) & is.numeric(above_merged_m)))) %>% 
  mutate_all(~replace_na(., 0)) %>%             # Remove Na
  decostand(method = "hellinger") %>%           # Relative abund
  colMeans() %>%                                # Mean relative abundances
  as.data.frame() %>% 
  mutate(rank = rank(-.)) %>% 
  setNames(c("ra", "rank"))




# Below Community Evenness ----
# Mustard Subsetting
below_merged_m <- subset(below_merged, plot_type == "mustard")

below_merged_relabun_m <- below_merged_m %>%  
  select(where(function(above_merged_m) any(!is.na(above_merged_m) & is.numeric(above_merged_m)))) %>% 
  mutate_all(~replace_na(., 0)) %>%             # Remove Na
  decostand(method = "hellinger") %>%           # Relative abund
  colMeans() %>%                                # Mean relative abundances
  as.data.frame() %>% 
  mutate(rank = rank(-.)) %>% 
  setNames(c("ra", "rank"))






# Fallow Subsetting
below_merged_f <- subset(below_merged, plot_type == "fallow")

below_merged_relabun_f <- below_merged_f %>%  
  select(where(function(above_merged_m) any(!is.na(above_merged_m) & is.numeric(above_merged_m)))) %>% 
  mutate_all(~replace_na(., 0)) %>%             # Remove Na
  decostand(method = "hellinger") %>%           # Relative abund
  colMeans() %>%                                # Mean relative abundances
  as.data.frame() %>% 
  mutate(rank = rank(-.)) %>% 
  setNames(c("ra", "rank"))

