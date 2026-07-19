# DBSLGM_code
# Code used for the manuscript on "Design-based spatial modelling of malaria incidence among children in Nigeria"
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
