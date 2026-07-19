# DBSLGM_code
# Code used for the manuscript on "Design-based spatial modelling of malaria incidence among children in Nigeria using Latent Gaussian model"
suppressMessages(library(INLA))
suppressMessages(library(haven))
suppressMessages(library(sf))
suppressMessages(library(tidyverse))
suppressMessages(library(tmap))
suppressMessages(library(excursions))
suppressMessages(library(spdep))
suppressMessages(library("rgdal"))
suppressMessages(library(mice))
suppressMessages(library(correlation))
suppressMessages(library(car))
suppressMessages(library(naniar))
suppressMessageslibrary(correlation))

# VARIABLE SELECTION AND DATA CLEANING
nmis2021 <- read_sav("NGKR81FL.sav")
nmis2021clus <- nmis2021 %>% select("V001","V002","ML501","V012","V106","V113","V459","V127","V128","V129","V171B","V218",
	"B4","B8","ML101","V190","V025","V005","H71","V024","V136","V137","V150","V151","V152","V160","V169A","V169C","V460","V461","V238","H22","ML503")
names(nmis2021clus)
names(nmis2021clus) <- c("EA_num","HH_num","messg","mth_age","mth_edu","sour_wtr","mos_net","flr_mt","wll_mt",
	"rf_mt","frq_int","lvg_chd","chd_sx","chd_age","typ_net","wt_ind","resd","wts","chd_mal","state","hh_no","no_chd5","hh_relat","sex_hh","age_hh","share_toilt","mb_phn","smrt_phn","chdslpmsqnt","slptmsqnt","birth_3yr","fevr_2wk","mal_prvtd")

#REMOVES ROWS WITH NAS in state and resd
nmis2021clus <-nmis2021clus %>% filter(!is.na(state))
nmis2021clus <-nmis2021clus %>% filter(!is.na(resd))

nmis2021clus[,c(3,5:11,13,15:17,19,20,23,24,26:30,32,33)] <- as_factor(nmis2021clus[,c(3,5:11,13,15:17,19,20,23,24,26:30,32,33)])

#Percentage of missing values
#Number of missing vaues
missingcount <- colSums(is.na(nmis2021clus))
#% of missing values
missingpercent <- colMeans(is.na(nmis2021clus))*100

#MAP MALARIA PREVALENCE
#Create IDs
nmis2021clus$area_id <- as.numeric(factor(nmis2021clus$state))
nmis2021clus$psu_id <- as.numeric(factor(nmis2021clus$EA_num))
nmis2021clus$strata_id <- as.numeric(factor(nmis2021clus$resd))
#SPATIAL exploration of malaria prevalence
Nig<- readOGR(".", "gadm36_NGA_1")

names(Nig)[names(Nig) == "NAME_1"] = "state" # rename column manually
names(Nig)
Nig@data$state[Nig@data$state == "Federal Capital Territory"] <- "FCT"
Nig@data$state[Nig@data$state == "Nassarawa"] <- "Nasarawa"

nb <- poly2nb(Nig)
#head(nb)
nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")

#MORAN'S I
#Aggregate the data to state level to count Yes(1) those with malaria
nmis2021clus4$state <- as_factor(nmis2021clus4$state)
state_data <- nmis2021clus4 %>% 
		 group_by(state)%>%
		summarise(chd_mal = sum(chd_mal == 1, na.rm=TRUE))
head(state_data)
weights <- nb2listw(nb, style="W")

#Morans I statistics
moran_result <- moran.test(state_data$chd_mal,weights)
print(moran_result)

#Spatial prevalence
#Merge spatial data
#state_data <- Nig$state
Nig.prev <- merge(Nig, state_data, by = "state", all.x = TRUE)
Nig.prev_sf <- st_as_sf(Nig.prev, cords =c("Longitude", "Latitude"), crs = 4326)
Nig.prev_sf <- st_transform(Nig.prev_sf,crs=4326)
tmap_mode("plot")

map1 <- tm_shape(Nig.prev_sf)+ tm_fill("chd_mal", style="pretty", pallete="viridis", title="Child Malaria Incidence")+
	tm_borders()+tm_layout(inner.margins=c(0.05,0.05,0.1,0.05),legend.text.size=0.75, legend.title.size=1.0, frame=FALSE)+ tm_graticules(lines=FALSE, labels.size=0.8)+ tm_compass(position=c("left","bottom"))+tm_scale_bar(text.size=0.75,position=c("centre","top"))+tm_text("state", size=0.7, xmod=0.5,ymod=0.5)
map1

#Multple Imputation
#Convert weight to decimals and Normalise
nmis2021clus <- nmis2021clus %>%
	mutate(wts_dec= wts/1000000)

nmis2021clus <-nmis2021clus %>% select(!wts)#Remove raw weights
#weight scaling
nmis2021clus$w_scaled <- nmis2021clus$wts_dec/mean(nmis2021clus$wts_dec)

#Remove labels in SPSS here in R with zap_labels()
nmis2021clus <- zap_labels(nmis2021clus)
meth <- make.method(nmis2021clus)
meth["smrt_phn"] <- "logreg"
meth["share_toilt"] <- "polyreg"
meth["chd_age"] <- "pmm"
meth["fevr_2wk"] <- "polyreg"
meth["chdslpmsqnt"] <- "polyreg" 
impp <- mice(nmis2021clus, m=10, method=meth, maxit=5, ridge = 0.01,seed=12345)

#Correlation
#Mixed-type association matrix for each imputed dataset
m <- impp$m
assoc.list <- vector("list", m)
for(i in 1:m){
dat <- complete(impp, i)
assoc.list[[i]] <- correlation(
dat[, c("mth_age","lvg_chd","chd_age","hh_no","no_chd5","age_hh","birth_3yr","messg","mth_edu","sour_wtr","mos_net", "flr_mt","wll_mt","rf_mt","frq_int","chd_sx","typ_net","wt_ind","resd","hh_relat","sex_hh","share_toilt", "mb_phn", "smrt_phn","chdslpmsqnt","slptmsqnt","fevr_2wk","mal_prvtd")])}

#pool assoc estimates
assoc.df <- bind_rows(lapply(assoc.list, as.data.frame),.id="Imputation")
assoc.summary <- assoc.df %>% group_by(Parameter1,Parameter2)%>%summarise(Mean=mean(r, na.rm=TRUE), SD=sd(r, na.rm=TRUE),Lower=quantile(r,0.025),Upper=quantile(r,0.975),.groups="drop")

#GVIF for each imputation
gvif.list <- vector("list",m)
for(i in 1:m){
dat <- complete(impp, i)
fit <- glm(chd_mal ~ mth_age +lvg_chd+chd_age+hh_no+no_chd5+age_hh+birth_3yr+messg+mth_edu+sour_wtr+mos_net+flr_mt+wll_mt+rf_mt+frq_int+chd_sx+typ_net+wt_ind+resd+hh_relat+sex_hh+share_toilt+mb_phn+smrt_phn+chdslpmsqnt+slptmsqnt+fevr_2wk+mal_prvtd, family=binomial, data=dat)
gvif.list[[i]] <- vif(fit)}

#Pool GVIF across imputations
gvif.df <- bind_rows(lapply(gvif.list, function(x){as.data.frame(x) |>tibble::rownames_to_column("Variable")}),
.id="Imputation")

#Summarise GVIF
gvif.summary<- gvif.df%>% group_by(Variable)%>% summarise(Mean_GVIF =mean(GVIF),SD=sd(GVIF), Minimum=min(GVIF),Maximum=max(GVIF),.groups="drop")

#Compute adjusted GVIF
gvif.summary <- gvif.df |> mutate(Adjusted_GVIF = GVIF^(1/(2*Df)))

#average adjusted GVIF
Adjusted.summary <- gvif.summary |> group_by(Variable) |> summarise(Mean_adjusted_GVIF = mean(Adjusted_GVIF), SD = sd(Adjusted_GVIF),.groups="drop")

#DBSLGM FITTING ON IMPUTED DATASETS
imp_list <- complete(imppp, action="all")

#ensure response is numeric
imp_list <- lapply(imp_list, function(dat){
dat$chd_mal <- ifelse(dat$chd_mal=="Yes",1,0)
dat})

#Create IDs
for(i in seq_along(imp_list)){
imp_list[[i]]$area_id <- as.numeric(factor(imp_list[[i]]$state))
imp_list[[i]]$psu_id <- as.numeric(factor(imp_list[[i]]$EA_num))
imp_list[[i]]$strata_id <- as.numeric(factor(imp_list[[i]]$resd))}

#specify PC priors
hyper.iid <- list(prec=list(prior="pc.prec", param=c(1,0.01)))
hyper.bym <- list(prec=list(prior="pc.prec", param=c(1,0.01)), phi=list(prior="pc",param=c(0.5,2/3)))

form_bym2 <- chd_mal ~ mth_age +lvg_chd+chd_age+hh_no+no_chd5+age_hh+messg+mth_edu+chd_sx+wt_ind+
 f(area_id, model="bym2", graph=g, scale.model=TRUE, hyper=hyper.bym)+
  f(psu_id, model="iid", hyper=hyper.iid)+
 f(strata_id, model="iid", hyper=hyper.iid)

#fit DBSLGM with BYM2 prior
models_bym2 <- vector("list", length(imp_list))
for(i in seq_along(imp_list)){
cat("BYM2 : Imputation", i, "\n")
models_bym2[[i]] <- inla(form_bym2, family="binomial", data=imp_list[[i]], Ntrials=1, weights=imp_list[[i]]$w_scaled, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE,config=TRUE)) }

#Average DIC
mean(sapply(models_bym2, function(x)x$dic$dic))

#Average WAIC
mean(sapply(models_bym2, function(x)x$waic$waic))

#Average LPML
LPML <- sapply(models_bym2, function(x)sum(log(x$cpo$cpo),na.rm=TRUE))
mean(LPML)

#ROC for each imputation
library(pROC)
roc.list <- lapply(1:2, function(i){
dat <- imp_list[[i]]
pred <- models_bym2[[i]]$summary.fitted.values$mean
roc(dat$chd_mal, pred)})
#Average AUC
mean(sapply(roc.list,auc))

#Posterior Predictive checks PIT
pit <- lapply(models_bym2,function(x)x$cpo$pit)
#PIT Histogram
pit.all <- unlist(lapply(models_bym2, function(x)x$cpo$pit))
pit.all <- pit.all[!is.na(pit.all)]
hist(pit.all, breaks=seq(0,1,by=0.05), probability=TRUE, col="lightblue", border="white",main="Pooled PIT Histogram Across Imputations",xlab="PIT")
abline(h=1, col="red",lwd=2,lty=2)

#Posterior Predictive Samples, I used 300 samples
samples <- lapply(models_bym2, function(x)inla.posterior.sample(1000,x))
#simluated replicated outcomes
sim.data <- lapply(1:2, function(i){
dat <- imp_list[[i]]
pred <- models_bym2[[i]]$summary.fitted.values$mean
Yrep <- matrix(rbinom(length(pred)*1000,1,rep(pred,1000)),nrow=length(pred))
Yrep})

#Posterior Predictive p-values
obs.stat <- lapply(1:2, function(i){
dat <- imp_list[[i]]
sum(dat$chd_mal)})
#replicated statistic
rep.stat <- lapply(sim.data,colSums)
#Bayesian posterior predictive p-value
ppp <- mapply(function(obs, rep){
mean(rep>=obs)},
obs.stat,rep.stat)
mean(ppp)

#Residual MORAN’S I for each imputation
lw <- nb2listw(nb, style="W")
moran.list <- lapply(1:2, function(i){
dat <- imp_list[[i]]
dat$residual <- dat$chd_mal-models_bym2[[i]]$summary.fitted.values$mean
state.res <- aggregate(residual~state, mean, data=dat)
#state$residual <- state.res$residual
moran.test(state.res$residual,lw)})

moran_residual <- data.frame(
	Imputation =1:2,
	Moran_I=sapply(moran.list, function(x)x$estimate[1]),
	Expected=sapply(moran.list, function(x)x$estimate[2]),
	Variance=sapply(moran.list, function(x)x$estimate[3]),
	Pvalue= sapply(moran.list, function(x)x$p.value))
moran_residual

#POoled Moran's I
summary_after <- data.frame(
	Mean_Moran =mean(moran_after$Moran_I),
	SD_Moran=sd(moran_after$Moran_I),
	Median_pv=median(moran_after$Pvalue))
summary_after

#CLassification Table
library(caret)
classification.list <- lapply(1:2, function(i){
dat <- imp_list[[i]]
pred <- models_bym2[[i]]$summary.fitted.values$mean
pred.class <- ifelse(pred >= 0.50,1,0)
confusionMatrix(factor(pred.class), factor(dat$chd_mal),positive="1")})

#Sensitivity and Specificty
sensitivity <- sapply(classification.list, function(x)x$byClass["Sensitivity"])
specificity <- sapply(classification.list, function(x)x$byClass["Specificity"])
mean(sensitivity)
mean(specificity)

#Posterior Predictive Probabilities
posterior.prob <- lapply(models_bym2, function(x)x$summary.fitted.values)
posterior.prob[[1]]$mean
posterior.prob[[1]][,c("0.025quant","0.975quant")]

#Posterior Predictive Residuals
residual.list <- lapply(1:2, function(i){
dat <- imp_list[[i]]
pred <- models_bym2[[i]]$summary.fitted.values$mean
dat$chd_mal-pred})

#Calibration Plot
pred_list <- lapply(models_bym2, function(mod){
	mod$summary.fitted.values$mean})
#pooled
pred_matrix <- do.call(cbind,pred_list)
pooled_pred <- rowMeans(pred_matrix)
#Pooled observation
obs <- imp_list[[1]]$chd_mal
#Calibration data
calibration <- data.frame(
	obs=obs, pred=pooled_pred)
calibration$decile <-ntile(calibration$pred,10)
#Observed and pred prob by decile
calib <- calibration %>% group_by(decile)%>%
	summarise(Observed=mean(obs),
	Predicted=mean(pred),n=n(),.groups="drop")

#plot calibration curve
ggplot(calibration, aes(x=pred, y=obs)) +
	geom_smooth(method="loess",se=TRUE,color="blue")+
	geom_abline(intercept=0,
			slope=1,
			linetype="dashed",
			colour="red")+
	coord_equal(xlim=c(0,1),ylim=c(0,1))+
	labs(title="LOESS Calibration Curve",
		x="Predicted Probabilty",
		y="Observed Malaria Incidence")+theme()

ggplot(calib, aes(x=Predicted, y=Observed)) +
	geom_point(size=3)+ geom_line()+
	geom_abline(intercept=0,
			slope=1,
			linetype="dashed",
			colour="red")+
	coord_equal(xlim=c(0,1),ylim=c(0,1))+
	labs(title="calibration Plot for Pooled Predictions",
		x="Mean Predicted Probabilty",
		y="Observed Proportion")+theme()

#Brier score
brier <- mean(obs - pooled_pred)^2
brier

#Extract fixed effects odd ratio
betas <- lapply(models_bym2, function(res){
	res$summary.fixed[,"mean"]|>exp()})
vars <- lapply(models_bym2, function(res){
	res$summary.fixed[,"sd"]*2})

#convert to matrix
beta_mat <- do.call(rbind, betas)
var_mat <- do.call(rbind, vars)

#Pool results for each parameter
m <- nrow(beta_mat)
#mean estimate
beta_bar <- colMeans(beta_mat)

#within imputation variance
W <- colMeans(var_mat)

#between imputation variance
B <- apply(beta_mat, 2, var)

#total variance
T <- W +(1 +1/m)*B

#standard error
se <- sqrt(T)

#confidence interval
lower <- beta_bar - 1.96*se
upper <- beta_bar + 1.96*se

#Final Output
results_table <- data.frame(Estimate = beta_bar, SE =se, lower95 = lower, upper95  =upper)

#Extract spatial random effects
spatial_means <- lapply(models_bym2, function(res){
	res$summary.random$area_id$mean
	})
spatial_vars <- lapply(models_bym2, function(res){
	res$summary.random$area_id$sd^2
	})
#convert to matrix spatial
spatial_mean_mat <- do.call(cbind, spatial_means)
spatial_vars_mat <- do.call(cbind, spatial_vars)
#Pool Spatial results for each parameter
m <- ncol(spatial_mean_mat)
#state mean across imputations
state_mean <- rowMeans(spatial_mean_mat)
#within imputation variance
W_state <- rowMeans(spatial_vars_mat)
#between imputation variance
B_state <- apply(spatial_vars_mat, 1, var)
#total variance
T_state <- W_state +(1 +1/m)*B_state
#standard error
se_state <- sqrt(T_state)
#confidence interval
state_lower <- state_mean - 1.96*se_state
state_upper <- state_mean + 1.96*se_state
#Final Output
#Spatial effects
state_results <- data.frame(sn = 1:length(state_mean), mean=state_mean, SE =se_state, lower = state_lower, upper  =state_upper)
state_results

#Compute relative risk per imputation
rr_list <- lapply(models_bym2, function(res){
	exp(res$summary.random$area_id$mean)})
#compute Exceedance probabilities
rr.threshold <- 1
eta.threshold <- log(rr.threshold)
prob_list <- lapply(models_bym2, function(mod){
	sapply(mod$marginals.random$area_id,
	function(m){
	1-inla.pmarginal(eta.threshold,m)})})
#average across imputation
rr_mat <- do.call(cbind, rr_list)
prob_mat <- do.call(cbind, prob_list)
rr_pooled <- rowMeans(rr_mat)
prob_pooled <- rowMeans(prob_mat)

#Mapping data
state_results$state <- Nig$state
state_results$RR <- rr_pooled
state_results$ExcProb <- prob_pooled
#Define clusters
state_results$clusters <- state_results$ExcProb > 0.90
str(state_results)

#Aggregate the data to state level 
state_data <- state_results %>% 
		 group_by(state)%>%
		summarise(randomEffect = mean(mean),
				RR_state = mean(RR), ExcProb_state=mean(ExcProb),
				Cluster=sum(clusters))

Nig.log <- merge(Nig, state_data, by = "state", all.x = TRUE)

Nig.log_sf <- st_as_sf(Nig.log, cords =c("Longitude", "Latitude"), crs = 4326)

#Merge spatial data
Nig.log_sf <- st_transform(Nig.log_sf,crs=4326)

#relative risk
map3 <- tm_shape(Nig.log_sf)+ tm_fill("RR_state", style="cont", pallete="viridis", title="Relative Risk")+
	tm_borders()+tm_layout(inner.margins=c(0.05,0.05,0.1,0.05),legend.text.size=0.75, legend.title.size=1.0, frame=FALSE)+ tm_graticules(lines=FALSE, labels.size=0.8)+ tm_compass(position=c("left","bottom"))+tm_scale_bar(text.size=0.75,position=c("centre","top"))+tm_text("state", size=0.7, xmod=0.5,ymod=0.5)
map3

#Map Exceedance Probabilites
map4 <- tm_shape(Nig.log_sf)+ tm_fill("ExcProb_state", style="cont", pallete="-RdYlBu", title="DBSLGM Exceedance Probability")+
	tm_borders()+tm_layout(inner.margins=c(0.05,0.05,0.1,0.05),legend.text.size=0.75, legend.title.size=1.0, frame=FALSE)+ tm_graticules(lines=FALSE, labels.size=0.8)+ tm_compass(position=c("left","bottom"))+tm_scale_bar(text.size=0.75,position=c("centre","top"))+tm_text("state", size=0.7, xmod=0.5,ymod=0.5)

#Map Bayesian Clusters
Nig.log_sf$Cluster <- factor(Nig.log_sf$Cluster,levels=c("1","0"),labels=c("Hotspot","Non-hotspot"))

map5 <- tm_shape(Nig.log_sf)+ tm_fill("Cluster", pallete=viridis::viridis(2,begin=0.2,end=0.8), title="Malaria Risk Clusters")+
	tm_borders()+tm_layout(inner.margins=c(0.05,0.05,0.1,0.05),legend.text.size=0.75, legend.title.size=1.0, frame=FALSE)+ tm_graticules(lines=FALSE, labels.size=0.8)+ tm_compass(position=c("left","bottom"))+tm_scale_bar(text.size=0.75,position=c("centre","top"))+tm_text("state", size=0.7, xmod=0.5,ymod=0.5)
map5

########################################################
## Bayesian Design-based Spatial Logistic Regression####
## Sensitivity Analyses################################# 

########################
####BYM, BYM2, Besag####
########################
form_bym <- chd_mal ~ mth_age +lvg_chd+chd_age+hh_no+no_chd5+age_hh+messg+mth_edu+chd_sx+wt_ind+
 f(area_id, model="bym", graph=g, scale.model=TRUE)+
  f(psu_id, model="iid", hyper=hyper.iid)+
 f(strata_id, model="iid", hyper=hyper.iid)

form_bym2 <- chd_mal ~ mth_age +lvg_chd+chd_age+hh_no+no_chd5+age_hh+messg+mth_edu+chd_sx+wt_ind+
 f(area_id, model="bym2", graph=g, scale.model=TRUE)+
  f(psu_id, model="iid", hyper=hyper.iid)+
 f(strata_id, model="iid", hyper=hyper.iid)

form_besag <- chd_mal ~ mth_age +lvg_chd+chd_age+hh_no+no_chd5+age_hh+messg+mth_edu+chd_sx+wt_ind+
 f(area_id, model="besag", graph=g, scale.model=TRUE)+
  f(psu_id, model="iid", hyper=hyper.iid)+
 f(strata_id, model="iid", hyper=hyper.iid)

#fit BYM
models_bym <- vector("list", length(imp_list))
for(i in seq_along(imp_list)){
cat("BYM : Imputation", i, "\n")
models_bym[[i]] <- inla(form_bym, family="binomial", data=imp_list[[i]], Ntrials=1, weights=imp_list[[i]]$w_scaled, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE,config=TRUE)) }

#fit BYM2
models_bym2 <- vector("list", length(imp_list))
for(i in seq_along(imp_list)){
cat("BYM2 : Imputation", i, "\n")
models_bym2[[i]] <- inla(form_bym2, family="binomial", data=imp_list[[i]], Ntrials=1, weights=imp_list[[i]]$w_scaled, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE,config=TRUE)) }

#fit Besag
models_besag <- vector("list", length(imp_list))
for(i in seq_along(imp_list)){
cat("Besag : Imputation", i, "\n")
models_besag[[i]] <- inla(form_besag, family="binomial", data=imp_list[[i]], Ntrials=1, weights=imp_list[[i]]$w_scaled, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE,config=TRUE)) }

#Compare DIC
dic_compare <- data.frame(Imputation=1:length(imp_list), BYM=sapply(models_bym, function(x)x$dic$dic), BYM2=sapply(models_bym2, function(x)x$dic$dic), Besag=sapply(models_besag, function(x)x$dic$dic))
dic_compare

#Compare waic
waic_compare <- data.frame(Imputation=1:length(imp_list), BYM=sapply(models_bym, function(x)x$waic$waic), BYM2=sapply(models_bym2, function(x)x$waic$waic), Besag=sapply(models_besag, function(x)x$waic$waic))
waic_compare

#Compare LPML
lpml <- function(model){
sum(log(model$cpo$cpo), na.rm=TRUE)}

lpml_compare <- data.frame(Imputation=1:length(imp_list), BYM=sapply(models_bym, lpml), BYM2=sapply(models_bym2, lpml), Besag=sapply(models_besag, lpml))
lpml_compare

#Average model performance
performance <- data.frame(Model=c("BYM","BYM2","Besag"), Mean_DIC=c(mean(dic_compare$BYM), mean(dic_compare$BYM2), mean(dic_compare$Besag)), Mean_WAIC=c(mean(waic_compare$BYM), mean(waic_compare$BYM2), mean(waic_compare$Besag)), Mean_LPML=c(mean(lpml_compare$BYM), mean(lpml_compare$BYM2), mean(lpml_compare$Besag)))
performance

#######################
####PC vs log-Gamma####
#######################
#PC prior
pc_prior <- list(prec=list(prior="pc.prec", param=c(1,0.01)))#=P(lambda>1)=0.01

#Alternative PC prior
pc_prior2 <- list(prec=list(prior="pc.prec", param=c(0.5,0.01)))#=P(lambda>0.5)=0.01

#Log-Gamma prior
loggamma_prior <- list(prec=list(prior="loggamma", param=c(1,0.0005)))

#Alternative Log-Gamma prior
loggamma_prior2 <- list(prec=list(prior="loggamma", param=c(0.5,0.005)))

#BYM2 with PC
form_pc <- chd_mal ~ mth_age +lvg_chd+chd_age+hh_no+no_chd5+age_hh+messg+mth_edu+chd_sx+wt_ind+
 f(area_id, model="bym2", graph=g, scale.model=TRUE, hyper=pc_prior)+
  f(psu_id, model="iid", hyper=hyper.iid)+
 f(strata_id, model="iid", hyper=hyper.iid)

#fit across all imputations
models_pc <- vector("list", length(imp_list))
for(i in seq_along(imp_list)){
models_pc[[i]] <- inla(form_pc, family="binomial", data=imp_list[[i]], Ntrials=1, weights=imp_list[[i]]$w_scaled, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE,config=TRUE)) }

#BYM2 with logGamma
form_lg <- chd_mal ~ mth_age +lvg_chd+chd_age+hh_no+no_chd5+age_hh+messg+mth_edu+chd_sx+wt_ind+
 f(area_id, model="bym2", graph=g, scale.model=TRUE, hyper= loggamma_prior)+
  f(psu_id, model="iid", hyper=hyper.iid)+
 f(strata_id, model="iid", hyper=hyper.iid)

#fit loggamma
models_lg <- vector("list", length(imp_list))
for(i in seq_along(imp_list)){
models_lg[[i]] <- inla(form_lg, family="binomial", data=imp_list[[i]], Ntrials=1, weights=imp_list[[i]]$w_scaled, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE,config=TRUE)) }

#Compare DIC
dic_prior <- data.frame(Imputation=1:length(imp_list), PC=sapply(models_pc, function(x)x$dic$dic), LogGamma=sapply(models_lg, function(x)x$dic$dic))
dic_prior

#Compare waic
waic_prior <- data.frame(Imputation=1:length(imp_list), PC=sapply(models_pc, function(x)x$waic$waic), LogGamma=sapply(models_lg, function(x)x$waic$waic))
waic_prior

#Compare LPML
lpml <- function(model){
sum(log(model$cpo$cpo), na.rm=TRUE)}

lpml_prior <- data.frame(Imputation=1:length(imp_list), PC=sapply(models_pc, lpml), LogGamma=sapply(models_lg, lpml))
lpml_prior

#Compare posterior ratio
extract_fixed <- function(model){
model$summary.fixed}

fixed_pc <- lapply(models_pc, extract_fixed)
fixed_lg <- lapply(models_lg, extract_fixed)
#compare posterior variance of spatial effects
spatial_var <- data.frame(Imputation=1:length(imp_list), PC=sapply(models_pc, function(x)x$summary.hyperpar$mean[1]), LogGamma=sapply(models_lg, function(x)x$summary.hyperpar$mean[1]))
spatial_var

##############################################
#####Queen vs Rook neighbourhood structure####
############################################## 
nb_queen <- poly2nb(Nig, queen=TRUE)
nb_rook <- poly2nb(Nig, queen=FALSE)

nb2INLA("queen.graph", nb_queen)
nb2INLA("rook.graph", nb_rook)

g_queen <- inla.read.graph(filename = "queen.graph")
g_rook <- inla.read.graph(filename = "rook.graph")

#BYM2 for queen
form_queen <- chd_mal ~ mth_age +lvg_chd+chd_age+hh_no+no_chd5+age_hh+messg+mth_edu+chd_sx+wt_ind+
 f(area_id, model="bym2", graph=g_queen, scale.model=TRUE, hyper= loggamma_prior)+
  f(psu_id, model="iid", hyper=hyper.iid)+
 f(strata_id, model="iid", hyper=hyper.iid)

#fit queen
models_queen <- vector("list", length(imp_list))
for(i in seq_along(imp_list)){
models_queen[[i]] <- inla(form_queen, family="binomial", data=imp_list[[i]], Ntrials=1, weights=imp_list[[i]]$w_scaled, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE,config=TRUE)) }

#BYM2 for rook
form_rook <- chd_mal ~ mth_age +lvg_chd+chd_age+hh_no+no_chd5+age_hh+messg+mth_edu+chd_sx+wt_ind+
 f(area_id, model="bym2", graph=g_rook, scale.model=TRUE, hyper= loggamma_prior)+
  f(psu_id, model="iid", hyper=hyper.iid)+
 f(strata_id, model="iid", hyper=hyper.iid)

#fit queen
models_rook <- vector("list", length(imp_list))
for(i in seq_along(imp_list)){
models_rook[[i]] <- inla(form_rook, family="binomial", data=imp_list[[i]], Ntrials=1, weights=imp_list[[i]]$w_scaled, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE,config=TRUE)) }

#Compare DIC
dic_nb <- data.frame(Imputation=1:length(imp_list), Queen=sapply(models_queen, function(x)x$dic$dic), Rook=sapply(models_rook, function(x)x$dic$dic))
dic_nb

#Compare waic
waic_nb <- data.frame(Imputation=1:length(imp_list), Queen=sapply(models_queen, function(x)x$waic$waic), Rook=sapply(models_rook, function(x)x$waic$waic))
waic_nb

#Compare LPML
lpml <- function(model){
sum(log(model$cpo$cpo), na.rm=TRUE)}

lpml_nb <- data.frame(Imputation=1:length(imp_list), Queen=sapply(models_queen, lpml), Rook=sapply(models_rook, lpml))
lpml_nb

#Compare posterior spatial risk
queen_rr <- lapply(models_queen, function(x)x$summary.random$state$mean)
rook_rr <- lapply(models_rook, function(x)x$summary.random$state$mean)

#compare exceedance prob
threshold <- 0

queen_exc <- lapply(models_queen, function(mod){
	sapply(mod$marginals.random$area_id,
		function(m){
		1-inla.pmarginal(threshold,m)})
	})

rook_exc <- lapply(models_rook, function(mod){
	sapply(mod$marginals.random$area_id,
		function(m){
		1-inla.pmarginal(threshold,m)})
	})

#Aggrement in hostpot classification of choice of 0.95
library(irr)
queen_hot <-as.integer(unlist(queen_exc[[1]]) > 0.95)
rook_hot <-as.integer(unlist(rook_exc[[1]]) > 0.95)

kappa2(data.frame(queen_hot, rook_hot))

#################################### 
#####Hyper-prior specifications#####
#################################### 
#default PC
pc1 <- list(prec=list(prior="pc.prec",param=c(1,0.01)))
#more informative
pc2 <- list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
#more difuse
pc3 <- list(prec=list(prior="pc.prec",param=c(2,0.01)))

#Log-Gamma 1
lg1 <- list(prec=list(prior="loggamma",param=c(1,0.0005)))
#Log-Gamma 1
lg2 <- list(prec=list(prior="loggamma",param=c(1,0.01)))
#Log-Gamma 1
lg3 <- list(prec=list(prior="loggamma",param=c(0.5,0.005)))

#Prior for BYM2 mixing parameter phi
phi1 <- list(phi=list(prior="pc",param=c(0.5,2/3)))
phi2 <- list(phi=list(prior="pc",param=c(0.5,0.5)))
phi3 <- list(phi=list(prior="pc",param=c(0.8,0.5)))

#list of all hyperpriors
hyper_list <- list(PC1 = list(prec=pc1$prec), PC2=list(prec=pc2$prec), PC3=list(prec=pc3$prec), LG1=list(prec=lg1$prec),LG2= list(prec=lg2$prec),LG3= list(prec=lg3$prec))

#Automatically fit every hyperprior
models_hyper <- list()
for(h in names(hyper_list)){
cat("Running",h, "\n")
models_hyper[[h]] <- vector("list",length(imp_list))
for(i in seq_along(imp_list)){
formula_hyper <- chd_mal ~ mth_age +lvg_chd+chd_age+hh_no+no_chd5+age_hh+messg+mth_edu+chd_sx+wt_ind+
 f(area_id, model="bym2", graph=g_rook, scale.model=TRUE, hyper= hyper_list[[h]])+
  f(psu_id, model="iid", hyper=hyper.iid)+
 f(strata_id, model="iid", hyper=hyper.iid)
models_hyper[[h]][[i]]<- inla(formula_hyper, family="binomial", data=imp_list[[i]], Ntrials=1, weights=imp_list[[i]]$w_scaled, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE,config=TRUE))}}

#Compare DIC
dic_hyper <- data.frame(Imputation=1:length(imp_list))
	for(h in names(models_hyper)){
dic_hyper[[h]] <- sapply(models_hyper[[h]], function(x)x$dic$dic)}
dic_hyper

#Compare waic
waic_hyper <- data.frame(Imputation=1:length(imp_list))
	for(h in names(models_hyper)){
waic_hyper[[h]] <- sapply(models_hyper[[h]], function(x)x$waic$waic)}
waic_hyper

#Compare LPML
lpml <- function(model){
sum(log(model$cpo$cpo), na.rm=TRUE)}

lpml_hyper <- data.frame(Imputation=1:length(imp_list))
for(h in names(models_hyper)){
lpml_hyper[[h]] <- sapply(models_hyper[[h]], lpml)}
lpml_hyper

#Average performance
performance <- data.frame(Prior=names(models_hyper), Mean_DIC=sapply(dic_hyper[-1],mean),
Mean_WAIC=sapply(waic_hyper[-1],mean),
Mean_LPML=sapply(lpml_hyper[-1],mean))
performance

#Posterior stability – Pooled posterior odds ratios
fixed_effects <- lapply(models_hyper, function(modlist){
lapply(modlist, function(x)x$summary.fixed)})

#Hyperparameter comparions
hyper_summary <- lapply(models_hyper, function(modlist){
lapply(modlist, function(x)x$summary.hyperpar)})

####################################### 
###Exceedance Probability Thresholds###
####################################### 
rr.threshold <- 1
eta.threshold <- log(rr.threshold)

#Extract exceedance prob
exc_prob <- lapply(models_bym2, function(mod){
	sapply(mod$marginals.random$area_id,
	function(m){
	1-inla.pmarginal(eta.threshold,m)})})

#Pool exceedance probability across imputations
pooled_exc <- Reduce("+", exc_prob)/length(exc_prob)
#posterior SD
exc_sd <- apply(do.call(cbind,exc_prob),1,sd)

#Define candidate thresholds
thresholds <- c(0.70, 0.80, 0.90, 0.95, 0.99)

#Classify hotspot
hotspots <- lapply(thresholds, function(th){
as.integer(pooled_exc >= th)})

#Assign names
names(hotspots) <- paste0("EP", thresholds)

#Number of hotspots
n.hot <- sapply(hotspots, sum)
n.hot

#Stability of classification
hot.matrix <- do.call(cbind, hotspots)#each row is state and col threshold

#Freq of being classified
stability <- rowMeans(hot.matrix)
stability

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




