# DBSLGM_code
# Code used for the manuscript on "Design-based spatial modelling of malaria incidence among children in Nigeria"

#Create table
table.threshold <- data.frame(Threshold =thresholds, Hotspots=n.hot)

#Identify robust hotspot regardless of threshold
robust <- which(stability ==1)

#Identify unstable hotspot, changing with threshold
unstable <- which(stability >0 & stability < 1)

#Pooled posterior relative risks
RR <- lapply(models_bym2, function(mod)exp(mod$summary.random$area_id$mean))

#Pooled RR
RR.pool <- Reduce("+", RR)/length(RR)

#Sensitivity figure
plot(thresholds, n.hot, type="b", pch=19, xlab="Exceedance Threshold", ylab="Number of Hotspots")

#Create maps
Nig$EP <- pooled_exc

Nig$Hot95 <- hotspots[[4]]




