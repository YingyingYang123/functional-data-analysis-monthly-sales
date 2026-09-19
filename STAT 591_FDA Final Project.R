#### STAT 591: FDA Final Project ####
## packages
library(tidyr)
library(dplyr)
library(fda)
library(lubridate) # Date format
library(zoo)
library(viridisLite) # adjust color

# Data
df = read.csv("C:/Users/Yingying Yang/Downloads/STAT591 FDA/STAT 591_FDA Final Project/train.csv",header=TRUE)
# Date format
df$Order.Date <- as.Date(df$Order.Date, format = "%d/%m/%Y")
df$Year <- year(df$Order.Date)
df$Month <- month(df$Order.Date)

#### Region over time
df_Region = df[c('Region', 'Sales', 'Year', 'Month')]

# Aggregate monthly sales per region
df_monthly <- df_Region %>%
  group_by(Region, Year, Month) %>%
  summarise(Monthly_Sales = sum(Sales), .groups = 'drop') %>%
  arrange(Region, Year, Month)

# Create a numeric time index for FDA
df_monthly$Time <- as.numeric(as.yearmon(paste(df_monthly$Year, df_monthly$Month), "%Y %m"))

# Split data by Region
regions <- unique(df_monthly$Region)
colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd")
names(colors) <- regions

# Scatter plot
plot(df_monthly$Time, df_monthly$Monthly_Sales, col = colors[df_monthly$Region], pch = 16,
     main = "Monthly Sales by Region",
     xlab = "Time", ylab = "Monthly Sales")
legend("topleft", legend = regions, col = colors, lwd = 2, pch = 16)

# Another plot for adding smoothed curves
plot(df_monthly$Time, df_monthly$Monthly_Sales, col = colors[df_monthly$Region], pch = 16,
     main = "Monthly Sales by Region (Smoothed Curves)",
     xlab = "Time", ylab = "Monthly Sales")
# Create basis
optimal_nbasis_list <- list()

for (reg in regions) {
  df_reg <- df_monthly %>% filter(Region == reg)
  
  if (nrow(df_reg) == 0 || any(is.na(df_reg$Time)) || any(is.na(df_reg$Monthly_Sales))) {
    next
  }
  
  time_range <- range(df_reg$Time)
  nbasis_values <- 5:15
  cv_errors <- numeric(length(nbasis_values))
  
  for (i in seq_along(nbasis_values)) {
    basis <- create.bspline.basis(rangeval = time_range, nbasis = nbasis_values[i])
    fdPar <- fdPar(basis, lambda = 0.8)
    result <- smooth.basis(df_reg$Time, df_reg$Monthly_Sales, fdPar)
    preds <- eval.fd(df_reg$Time, result$fd)
    cv_errors[i] <- mean((df_reg$Monthly_Sales - preds)^2)
  }
  
  best_nbasis <- nbasis_values[which.min(cv_errors)[1]] 
  optimal_nbasis_list[[reg]] <- best_nbasis            
  
  best_basis <- create.bspline.basis(rangeval = time_range, nbasis = best_nbasis)
  fdobj <- smooth.basis(df_reg$Time, df_reg$Monthly_Sales, fdPar(best_basis, lambda = 1e-2))$fd
  
  lines(fdobj, col = colors[reg], lwd = 2)
}
legend("topleft", legend = regions, col = colors, lwd = 2, pch = 16)
optimal_nbasis_list


#### State over time
df_State <- df[c('State', 'Sales', 'Year', 'Month')]

# Aggregate monthly sales per state
df_monthly <- df_State %>%
  group_by(State, Year, Month) %>%
  summarise(Monthly_Sales = sum(Sales), .groups = 'drop') %>%
  arrange(State, Year, Month)

# Create a continuous time variable
df_monthly$Time <- as.numeric(as.yearmon(paste(df_monthly$Year, df_monthly$Month), "%Y %m"))

# Filter states with at least 5 months of data
state_counts <- df_monthly %>%
  group_by(State) %>%
  summarise(n = n()) %>%
  filter(n >= 5)

valid_states <- state_counts$State
states <- unique(valid_states)

# Assign colors using viridis
colors <- setNames(viridis(length(states), option = "A"), states)

# Scatter plot of monthly sales by state
plot(df_monthly$Time, df_monthly$Monthly_Sales,
     col = colors[df_monthly$State],
     pch = 16,
     main = "Monthly Sales by State",
     xlab = "Time", ylab = "Monthly Sales")
legend("topleft", legend = states, col = colors, lwd = 2, cex = 0.35, bty = "n")

# plot: smoothed curves by state
plot(df_monthly$Time, df_monthly$Monthly_Sales,
     col = colors[df_monthly$State],
     pch = 16,
     main = "Monthly Sales by State (Smoothed Curves)",
     xlab = "Time", ylab = "Monthly Sales")

optimal_nbasis_list <- list()

for (sta in states) {
  df_sta <- df_monthly %>% filter(State == sta)
  
  if (nrow(df_sta) == 0 || any(is.na(df_sta$Time)) || any(is.na(df_sta$Monthly_Sales))) {
    next
  }
  
  time_range <- range(df_sta$Time)
  nbasis_values <- 5:15
  cv_errors <- numeric(length(nbasis_values))
  
  for (i in seq_along(nbasis_values)) {
    basis <- create.bspline.basis(rangeval = time_range, nbasis = nbasis_values[i])
    fdPar <- fdPar(basis, lambda = 0.8)
    result <- smooth.basis(df_sta$Time, df_sta$Monthly_Sales, fdPar)
    preds <- eval.fd(df_sta$Time, result$fd)
    cv_errors[i] <- mean((df_sta$Monthly_Sales - preds)^2)
  }
  
  best_nbasis <- nbasis_values[which.min(cv_errors)[1]] 
  optimal_nbasis_list[[sta]] <- best_nbasis            
  
  best_basis <- create.bspline.basis(rangeval = time_range, nbasis = best_nbasis)
  fdobj <- smooth.basis(df_sta$Time, df_sta$Monthly_Sales, fdPar(best_basis, lambda = 1e-2))$fd
  
  lines(fdobj, col = colors[sta], lwd = 2)
}

legend("topleft", legend = states, col = colors, lwd = 2, cex = 0.35, bty = "n")

optimal_nbasis_list


#### functional format
df_wide <- df_monthly %>%
  mutate(Time = paste0(Year, "_", sprintf("%02d", Month))) %>%
  dplyr::select(State, Time, Monthly_Sales) %>%
  pivot_wider(names_from = Time, values_from = Monthly_Sales)

time_points <- sort(unique(df_monthly$Time))
state_names <- df_wide$State
sales_mat <- as.matrix(df_wide[,-1])  # remove State column
sales_mat[is.na(sales_mat)] <- 0

basis <- create.bspline.basis(rangeval = c(1, length(time_points)), nbasis = 10)
fd_obj <- Data2fd(argvals = 1:length(time_points), y = t(sales_mat), basisobj = basis)

fpca <- pca.fd(fd_obj, nharm = 10)  

prop_var <- fpca$varprop
cum_var <- cumsum(prop_var)

plot(cum_var, type = 'b', pch = 19,
     xlab = "Number of Principal Components",
     ylab = "Cumulative Proportion of Variance Explained",
     main = "Cumulative Variance Explained by FPCA")

K_selected <- which(cum_var >= 0.95)[1]

# K_selected fPC basis
par(mfrow = c(2, 2)) 
for (k in 1:K_selected) {
  plot(fpca$harmonics[k], 
       xlab = "Time", 
       ylab = paste0("fPC", k), 
       main = paste("Functional PC", k))
}

par(mfrow = c(1, 1))


set.seed(123)
# Elbow plot to choose optimal k
wss <- numeric(10)
for (k in 1:10) {
  kmeans_res <- kmeans(fpca$scores[, 1:4], centers = k, nstart = 25)
  wss[k] <- kmeans_res$tot.withinss
}

# Plot for Elbow Method
plot(1:10, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters K",
     ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Method for Choosing K")

# Choose optimal k
optimal_k <- 3
kmeans_final <- kmeans(fpca$scores[, 1:4], centers = optimal_k, nstart = 25)

# Cluster plot
plot(fpca$scores[, 1:2], col = kmeans_final$cluster, pch = 19,
     xlab = "PC1", ylab = "PC2", main = paste("Clustering Based on FPCA (k =", optimal_k, ")"))

cluster_df <- data.frame(
  State = state_names,
  Cluster = kmeans_final$cluster
)

# states in each cluster
split(cluster_df$State, cluster_df$Cluster)




