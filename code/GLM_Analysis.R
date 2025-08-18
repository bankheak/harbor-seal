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
require(ggplot2)
require(performance) # for overdispersion
require(AICcmodavg) #for AICc
require(glmmTMB) # for glmm
require(MASS) # for glm
require(DHARMa) # auto cor and zero-inflation
require(HH) # VIF


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
new.w.data$j.date<- as.integer(new.w.data$j.date)
new.w.data$j.date<- ifelse(new.w.data$j.date == min(new.w.data$j.date),0, new.w.data$j.date-min(new.w.data$j.date))
m.data$j.date<- as.integer(m.data$j.date)
m.data$j.date<- ifelse(m.data$j.date == min(m.data$j.date),0, m.data$j.date-min(m.data$j.date))
new.w.data$year <-as.numeric(format(new.w.data$date, format = "%Y"))
m.data$year <-as.numeric(format(m.data$date, format = "%Y"))


# Add sample ID
new.w.data$id<- seq_along(new.w.data[,1])
m.data$id<- seq_along(m.data[,1])

# Merge data
full.data<-merge(new.w.data,m.data,all = T)
full.data$site<- as.numeric(full.data$site=="waterfront")
summary(full.data)

## Organize data by id
full.data<- full.data[order(full.data$id), ]
full.data$obs_id<- seq_along(full.data[,1])
full.data<- full.data[order(full.data$obs_id), ]

## Create csv
write.csv(full.data,"full.data.csv")
full.data<- read.csv("full.data.csv")

# Check distribution
ggplot(full.data, aes(x=seals)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density()+stat_density(alpha=.2,adjust = 1, fill="#FF6666")+xlab("Number of Seals Hauled-out")+ylab("Density")+theme(panel.background = element_blank())

## Does the mean equal the variance?
var(full.data$seals)
mean(full.data$seals)
### Variance is way higher than the mean

## Negative binomial or poisson would be the best fit

# Check Temporal Autocorrelation
## Run acf for both direct and indirect effects
acf(full.data$seals[full.data$site == "waterfront"]) # Highly autocorrelated at lag of 1
acf(full.data$seals[full.data$site == "marina"]) # Autocorrelated at lag of 6

### Check Residuals of Month
fit1<- lm(seals~ month + year, data = full.data)
plot(full.data$month, fit1$res, pch=20, col="blue")
abline(h=0) # Add the horizontal line at 0

# Correct for autocorrelation with ar1()
fittedModel<- glmmTMB(seals ~ noise*site + tide + time +
                          ar1(as.factor(month) + 0 | year), data = full.data, 
        family = nbinom2, zi = ~ 1)
res = simulateResiduals(fittedModel)
testTemporalAutocorrelation(res, time = unique(full.data$obs_id))

# Test for homogeneity of variance
testCategorical(res, full.data$site)

## Check if overdispersion is detected with full model
check_overdispersion(fittedModel) 

### Check if data is zero inflated 
testZeroInflation(res) 
### Data is zero-inflated

# Run pairwise cor between all independent variables
## Cut-off is +/- 0.7
cor.matrix<-cor(full.data[,c(2:4,6:7,9)]) 
# Keep only large correlations in the same model
cor.matrix[abs(cor.matrix)< 0.7]<-NA
cor.matrix
## No colilinearity

###########################################################################
# PART 3: Run glmmTMB and AICc---------------------------------------------
  

fittedModel<- glmmTMB(count ~ avgnoise + Vmag + time + ar1(0 | month), 
                      data = full_data, family = poisson, zi = ~ 1) 

# Run glmmTMB
model0<- glmmTMB(seals ~ 1, data = full.data, 
                 family = nbinom2, zi = ~ 1)
summary(model0)
model1<- glmmTMB(seals ~ noise*site + tide + time +
                        ar1(as.factor(month) + 0 | year), data = full.data, 
                      family = nbinom2, zi = ~ 1)
summary(model1) 
model2<- glmmTMB(seals ~ noise*site + time +
                   ar1(as.factor(month) + 0 | year), data = full.data, 
                 family = nbinom2, zi = ~ 1)
summary(model2) 
model3<- glmmTMB(seals ~ noise*site + 
                   ar1(as.factor(month) + 0 | year), data = full.data, 
                 family = nbinom2, zi = ~ 1)
summary(model4) 
model4<- glmmTMB(seals ~ noise+site + 
                   ar1(as.factor(month) + 0 | year), data = full.data, 
                 family = nbinom2, zi = ~ 1)
summary(model4)

# Create list
models<-list(model0, model1, model2, model3, model4)

## Calculate AICc with glm of models
Modnames<- c("seals ~ 1",
             "seals ~ site*noise + tide + time",
             "seals ~ site*noise + time",
             "seals ~ site*noise",
             "seals ~ site+noise")

aictable<- aictab(models, modnames = Modnames,
                  second.ord = TRUE, 
                  nobs = NULL, c.hat = 1) # Looks like site*noise are the best predictors

# Find significant predictors and summary of model3
summary(model3) # noise between sites are significantly different


###########################################################################
# PART 4: Run Predict Function---------------------------------------------

# Use model to predict the response variable 
newdata <- data.frame(noise = mean(full.data$noise), 
                       site = factor(1:2, levels = 1:2, 
                                      labels = unique(full.data$site)))
m3<- glm.nb(seals ~ site * noise, data = full.data)
newdata$phat <- predict.glm(m3, newdata, type = "response")
newdata

# Interaction plot
interaction.plot(x.factor = full.data$noise, #x-axis variable
                 trace.factor = full.data$site, #variable for lines
                 response = full.data$seals, #y-axis variable
                 fun = median, #metric to plot
                 ylab = "Number of Seals Hauled-out",
                 xlab = "Noise Level (dB)",
                 col = c("pink", "blue"),
                 lty = 1, #line type
                 lwd = 2, #line width
                 trace.label = "Site")

# Plot 
newdata2 <- data.frame(
  noise = rep(seq(from = min(full.data$noise), to = max(full.data$noise), length.out = 100), 2),
  site = factor(rep(1:2, each = 100), levels = 1:2, labels =
                  unique(full.data$site)))

newdata2 <- cbind(newdata2, predict(m3, newdata2, type = "link", se.fit=TRUE))
newdata2 <- within(newdata2, {
  seals <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

write.csv(newdata2, "newdata2.csv")

ggplot(newdata2, aes(noise, seals)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = site), alpha = .25) +
  geom_line(aes(colour = site, linetype = site), linewidth = 1) +
  labs(x = "Noise Level (dB)", y = "Predicted Number of Seals Hauled-out") +
  theme(panel.background = element_blank())
