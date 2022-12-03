# Anthropogenic noise disturbance on harbor seals

###########################################################################
# Waterfront-Marina GLM analysis
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/habor-seal/data")

# Retrieve data
m.data<-read.csv("m.data.csv")
w.data<-read.csv("w.data.csv")

# Load packages
require(lme4)
require(ggplot2)
require(performance) #for check_collinearity()
require(MuMIn) # For dredge function
require(MASS) #for glm.nb()

###########################################################################
# PART 1: Check Location Significance at Waterfront ---------------------------------------------

head(w.data)

# Check t-test assumptions
## Is noise normal?
hist(w.data$noise[w.data$location == 1])
hist(w.data$noise[w.data$location == 2])
### Both are normal

## Equal variance?
var(w.data$noise[w.data$location == 1])
var(w.data$noise[w.data$location == 2])
### Equal variance

# Run t-test
t.test(w.data$noise[w.data$location == 1], w.data$noise[w.data$location == 2])
## Not significantly different
## Move on-to merging data sets

###########################################################################
# PART 2: Combine Waterfront and Marina Models and find distribution---------------------------------------------

# Combine location avg noise and total seals within dates
new.w.data<-read.csv("new.w.data.csv")

# Fix date
new.w.data$date<-as.Date(as.character(new.w.data$date),format = "%m/%d/%Y")
m.data$date<-as.Date(as.character(m.data$date),format = "%m/%d/%Y")

# Merge data
full.data<-merge(new.w.data,m.data,all = T)
summary(full.data)
full.data$time<-as.numeric(full.data$time)
full.data$month<-as.numeric(full.data$month)

# Run pairwise cor between all independent variables
## Cut-off is +/- 0.7
cor.matrix<-cor(full.data[,c(2:4,6:7)]) 
# Keep only large correlations in the same model
cor.matrix[abs(cor.matrix)< 0.7]<-NA
cor.matrix
## No colilinearity

# Check distribution
ggplot(full.data, aes(x=seals)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density()+stat_density(alpha=.2,adjust = 1, fill="#FF6666")+xlab("Number of Seals Hauled-out")+ylab("Density")+theme(panel.background = element_blank())
## Negative binomial would be the best fit

###########################################################################
# PART 3: Run GLM and AICs---------------------------------------------

# Create presence column
full.data$Presence<- ifelse(full.data$seals == 0, 0, 1)

# Change na.action
options(na.action = "na.fail")

# Fit model with all uncorrelated parameters
summary(all.parms<-glm(Presence~scale(j.date) + tide + time + site + noise + month, data = full.data, family = binomial))

# The dredge function fits all combinations
# of the variables in the all.parms model fit above
results<-dredge(all.parms)

# Grab best model based on lowest AICc
subset(results, delta == 0)
## Looks like month, noise, and site

# Check GLM
model<- glm.nb(seals ~ site*noise + month, data = full.data)
summary(model) #everything is significant


###########################################################################
# PART 4: Run Diagnostics---------------------------------------------

# Plot residuals by predicted values
plot(model$resid~model$fitted) # Plot residual vs. predicted values
plot(model$resid~full.data$noise)
plot(model$resid~full.data$month)
## All values are zero inflated
