# libraries ---------------------------------------------------------------
library(tidyverse)


# data --------------------------------------------------------------------
# Taken from Engemann, K., Sandel, B., Boyle, B., Enquist, B. J., Jørgensen, P. M., Kattge, J., McGill, B. J., Morueta-Holme, N., Peet, R. K., Spencer, N. J., Violle, C., Wiser, S. K. and Svenning, J.-C. (2016), 
# A plant growth form dataset for the New World. Ecology, 97: 3243. doi:10.1002/ecy.1569
GrowForm<-read.table("./data/base/GrowthForm_Final.txt",header=TRUE)

# traits data
# maximum plant height (m)
# SLA (cm2/g), 
# seed mass (mg)
# Leaf phosphorus and leaf nitrogen concentration per mass (Leaf N and Leaf P) (mg/g)
# wood density (mg/cm3).
Traits_BIEN<-read.csv("./data/base/2018_02_07_BIEN_trait_data.csv")


# Filter traits -----------------------------------------------------------
## Include only the six main trait levels
fun_traits<-c("whole plant leaf area per whole plant leaf dry mass",
              "seed mass","whole plant height","stem wood density",
              "leaf nitrogen content per leaf dry mass","leaf phosphorus content per leaf dry mass")

Traits_BIEN_sub<-
  Traits_BIEN %>%
  filter(trait_name%in%fun_traits)

Traits_BIEN_sub<-droplevels(Traits_BIEN_sub)

# Filter non numerical trait values -----------------------------------------
## Trait values such as "." or "*" and even "0"
Traits_BIEN_sub<-
  Traits_BIEN_sub %>% 
  filter(trait_value!="."&trait_value!="*"&trait_value!="0")

Traits_BIEN_sub<-droplevels(Traits_BIEN_sub)

## Exploring the number of traits with values per species
species_coverage<-
  Traits_BIEN_sub %>%
  group_by(scrubbed_species_binomial) %>% 
  summarise(N_traits=n_distinct(trait_name))

# Reshaping data frame ----------------------------------------------------
#Renaming trait factors
Traits_BIEN_sub$trait_name<-as.character(Traits_BIEN_sub$trait_name)
Traits_BIEN_sub$trait_name[which(Traits_BIEN_sub$trait_name=="stem wood density")]<-"Wood_density"
Traits_BIEN_sub$trait_name[which(Traits_BIEN_sub$trait_name=="leaf nitrogen content per leaf dry mass")]<-"Leaf_N"
Traits_BIEN_sub$trait_name[which(Traits_BIEN_sub$trait_name=="seed mass")]<-"Seed_mass"
Traits_BIEN_sub$trait_name[which(Traits_BIEN_sub$trait_name=="whole plant height")]<-"Height"
Traits_BIEN_sub$trait_name[which(Traits_BIEN_sub$trait_name=="leaf phosphorus content per leaf dry mass")]<-"SLA"
Traits_BIEN_sub$trait_name<-as.factor(Traits_BIEN_sub$trait_name)


# Number of observations per trait values in each species
# For example, Abarema jupunba has 16 values of stem_wood_density
Traits_BIEN_sub %>%
  group_by(scrubbed_species_binomial,trait_name) %>% 
  tally()

# Same but using the count function
Traits_BIEN_sub %>%
  count(scrubbed_species_binomial,trait_name)


Traits_BIEN_sub$trait_value_NU<-as.numeric(as.character(Traits_BIEN_sub$trait_value))
## Nas are produced in some fields where the trait value is characters
# such as "dead" "sacrificed"

# Calculate main trait values per species ---------------------------------
# Only one trait value per species
Mean_Traits_BIEN <-
  Traits_BIEN_sub %>% 
  group_by(scrubbed_species_binomial,trait_name) %>% 
  summarise(trait_value=mean(trait_value_NU,na.rm=TRUE))

var_names<-as.character(unique(Mean_Traits_BIEN$trait_name))

# Create trait values table, I use sum function to get rid of the NAs
# using mean also work since there is only one trait value per trait name
trait_df_wide<-
  Mean_Traits_BIEN  %>% 
  spread(trait_name,trait_value) %>% 
  group_by(scrubbed_species_binomial) %>%
  summarise_at(var_names,sum,na.rm=TRUE)

# Changing zeros to NAS
trait_df_wide<-as.data.frame(trait_df_wide)
trait_df_wide[trait_df_wide == 0]<-NA


# Include growth form -----------------------------------------------------
GrowForm_tmp<-
  GrowForm %>% 
  dplyr::select(FAMILY_STD,SPECIES_STD,GROWTHFORM_STD,GROWTHFORM_DIV)

Trait_BIEN_df<-merge(x=trait_df_wide,y=GrowForm_tmp,
                     by.x="scrubbed_species_binomial",
                     by.y="SPECIES_STD",
                     all.x=TRUE)

# Include a new Growth form with a more general classification
woody<-c("Tree","Liana","Shrub","Woody epiphyte")
Trait_BIEN_df$GROWTHFORM_GEN<-ifelse(Trait_BIEN_df$GROWTHFORM_STD%in%woody,"woody","herbaceous")
Trait_BIEN_df$GROWTHFORM_GEN[which(is.na(Trait_BIEN_df$GROWTHFORM_STD))]<-NA

Trait_BIEN_df %>% 
  group_by(GROWTHFORM_GEN) %>% 
  summarise(sp_number=n_distinct(scrubbed_species_binomial))

# Write clean datasets ----------------------------------------------------
write.csv(Trait_BIEN_df,"./data/processed/BIEN_trait_GrowthForm.csv")
