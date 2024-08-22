setwd("~/Desktop/r-summerproject/glucose-analysis")
library(dplyr)
library(lme4)
library(readr)
library(cowplot)
library(stringr)
library(car)
library(tidymodels)
library(boot)
library(readxl)
library

gluc.net <- read_excel("gluc-net-data-14aug.xlsx", 
                            sheet = "Sheet1")
alpha_diversity <- read_csv("~/Desktop/r-summerproject/glucose-analysis/a-div-physeq.csv")
alpha_diversity <- alpha_diversity %>%
  rename(sample.id = ...1)
alpha_diversity <- alpha_diversity %>%
  mutate(temperature = as.numeric(str_extract(sample.id, "\\d+")))


# 1. Exploratory scatterplots of experimental dataset: Chao1, Shannon, and Richness
chao.div.plot <- ggplot(alpha_diversity, aes(x = temperature, y = Chao1, group = temperature)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 6, color = "red", fill = "red") +
  theme_cowplot(12) +
  geom_jitter(shape = 16, position = position_jitter(0.2))

shannon.div.plot <- ggplot(alpha_diversity, aes(x = temperature, y = Shannon, group = temperature)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 6, color = "red", fill = "red") +
  theme_cowplot(12) +
  geom_jitter(shape = 16, position = position_jitter(0.2))

richness.div.plot <- ggplot(alpha_diversity, aes(x = temperature, y = Observed, group = temperature)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 6, color = "red", fill = "red") +
  theme_cowplot(12) +
  geom_jitter(shape = 16, position = position_jitter(0.2))


# 2. Simple linear regressions looking at diversity x temperature
lm.chao1 <- lm(formula = Chao1 ~ temperature, data = alpha_diversity)
summary(lm.chao1)
shapiro.test(alpha_diversity$Chao1) # Normality test: non-significant (normal)
ncvTest(lm.chao1) # Homoscedasticity: non-significant (not heteroscedastic)

lm.shannon <- lm(formula = Shannon ~ temperature, data = alpha_diversity)
summary(lm.shannon)
shapiro.test(alpha_diversity$Shannon) # Normality test: non-significant (normal)
ncvTest(lm.shannon) # Homoscedasticity: non-significant (not heteroscedastic)

lm.observed <- lm(formula = Observed ~ temperature, data = alpha_diversity)
summary(lm.observed)
shapiro.test(alpha_diversity$Observed) # Normality test: non-significant (normal)
ncvTest(lm.observed) # Homoscedasticity: non-significant (not heteroscedastic)


# 3. Bootstrapping Shannon 
shannon.lm.plot <- ggplot(alpha_diversity, aes(temperature, Shannon)) +
  geom_point() +
  geom_line(aes(y = predict(lm.shannon)))

set.seed(27) # get 10000 bootstrapped replicates
boots.shannon <- bootstraps(alpha_diversity, times = 10000, apparent = TRUE)
boots.shannon

fit.shannon.lm.boots <- function(split) { # create a function to run the linear model on bootstrapped replicates
  lm(formula = Shannon ~ temperature, analysis(split), start = list(k = 1, b = 0))
}

shannon.lm.boots <- boots.shannon %>%
  mutate(model = map(splits, fit.shannon.lm.boots),
         coef_info = map(model, tidy))

shannon.lm.boots.coefs <- shannon.lm.boots %>% # obtain coefficients from each bootstrapped model. 
  unnest(coef_info)

shannon.lm.perc.int <- int_pctl(shannon.lm.boots, coef_info) # obtain confidence intervals from the bootstrapping
shannon.lm.perc.int

shannon.lm.boots.aug <- shannon.lm.boots %>% # prepare bootstrapped data for plotting
  sample_n(200) %>%
  mutate(augmented = map(model, augment)) %>%
  unnest(augmented)

shannon.ci <- shannon.lm.boots.aug %>%
  group_by(temperature) %>%
  summarize(ci_lower = quantile(.fitted, 0.025),
            ci_upper = quantile(.fitted, 0.975))

# Plot with ggplot2
shannon.lm.boots.plot <- ggplot() + 
  geom_point(data = shannon.lm.boots.aug, aes(x = temperature, y = Shannon)) +
  geom_ribbon(data = shannon.ci, aes(x = temperature, ymin = ci_lower, ymax = ci_upper), 
              fill = "grey", alpha = 0.5) +
  geom_line(data = alpha_diversity, aes(x = temperature, y = predict(lm.shannon)), colour = "#DC4D01", linewidth = 1.5) +
  labs(x = "Temperature (°C)", y = "Shannon") +
  theme_cowplot(12)

# Print the plot
print(shannon.lm.boots.plot)

# 4. Bootstrapping Observed
observed.lm.plot <- ggplot(alpha_diversity, aes(temperature, Observed)) +
  geom_point() +
  geom_line(aes(y = predict(lm.observed)))

set.seed(27) # get 10000 bootstrapped replicates
boots.observed <- bootstraps(alpha_diversity, times = 10000, apparent = TRUE)
boots.observed

fit.observed.lm.boots <- function(split) { # create a function to run the linear model on bootstrapped replicates
  lm(formula = Observed ~ temperature, analysis(split), start = list(k = 1, b = 0))
}

observed.lm.boots <- boots.observed %>%
  mutate(model = map(splits, fit.observed.lm.boots),
         coef_info = map(model, tidy))

observed.lm.boots.coefs <- observed.lm.boots %>% # obtain coefficients from each bootstrapped model. 
  unnest(coef_info)

observed.lm.perc.int <- int_pctl(observed.lm.boots, coef_info) # obtain confidence intervals from the bootstrapping
observed.lm.perc.int

observed.lm.boots.aug <- observed.lm.boots %>% # prepare bootstrapped data for plotting
  sample_n(200) %>%
  mutate(augmented = map(model, augment)) %>%
  unnest(augmented)

observed.ci <- observed.lm.boots.aug %>%
  group_by(temperature) %>%
  summarize(ci_lower = quantile(.fitted, 0.025),
            ci_upper = quantile(.fitted, 0.975))

# Plot with ggplot2
observed.lm.boots.plot <- ggplot() + 
  geom_point(data = observed.lm.boots.aug, aes(x = temperature, y = Observed)) +
  geom_ribbon(data = observed.ci, aes(x = temperature, ymin = ci_lower, ymax = ci_upper), 
              fill = "grey", alpha = 0.5) +
  geom_line(data = alpha_diversity, aes(x = temperature, y = predict(lm.observed)), colour = "#DC4D01", linewidth = 1.5) +
  labs(x = "Temperature (°C)", y = "Observed OTUs") +
  theme_cowplot(12)

# Print the plot
print(observed.lm.boots.plot)
  
# 4. Plotting results from network analysis
gluc.net <- gluc.net %>%
  mutate(ave.path.len.diff = "YES") %>%
  mutate(nat.conn.diff = if_else(temperature > 20, "YES", "NO"))

# Get this data from the 100 random networks generated in the "gluc-net-met.R" file
b15 <- read_csv("~/Desktop/r-summerproject/glucose-analysis/globprop-15.csv")
b20 <- read_csv("~/Desktop/r-summerproject/glucose-analysis/globprop-20.csv")
b25 <- read_csv("~/Desktop/r-summerproject/glucose-analysis/globprop-25.csv")
b30 <- read_csv("~/Desktop/r-summerproject/glucose-analysis/globprop-30.csv")
combined <- rbind(b15, b20, b25, b30)
combined <- combined %>%
  rename(temperature = condition) %>%
  rename(nat.conn = natConnect1) %>%
  rename(ave.path.len = avPath1)

nc.boot.plot <- ggplot(combined, aes(x = temperature, y = nat.conn, group = temperature)) +
  geom_jitter(width = 1, color = "grey", alpha = 0.5)
nc.boot.plot <- nc.boot.plot +
  geom_point(data = gluc.net, aes(x = temperature, y = nat.conn), color = "#DC4D01", size = 3) +
  labs(x = "Temperature (°C)", y = "Natural connectivity") +
  theme_cowplot(12)
nc.boot.plot

ap.boot.plot <- ggplot(combined, aes(x = temperature, y = ave.path.len, group = temperature)) +
  geom_jitter(width = 1, color = "grey", alpha = 0.5)
ap.boot.plot <- ap.boot.plot +
  geom_point(data = gluc.net, aes(x = temperature, y = ave.path.len), color = "#DC4D01", size = 3) +
  labs(x = "Temperature (°C)", y = "Average path length") +
  theme_cowplot(12)
ap.boot.plot


