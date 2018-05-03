

# functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")

# data --------------------------------------------------------------------
# Presence data (included biomes)
load("./outputs/spPresence_biomes_all.RData")

# Growth form
Growth_form<-read.table("./data/base/GrowthForm_Final.txt",header = TRUE)

Growth_form$SPECIES_STD<-gsub(" ","_",Growth_form$SPECIES_STD)


spPresence_biomes<-spPresence %>% 
  select(Species,biomes) %>% 
  distinct

indx<-match(spPresence_biomes$Species, Growth_form$SPECIES_STD)
spPresence_biomes$GROWTHFORM_STD<-Growth_form$GROWTHFORM_STD[indx]


growth_biomes<-
  spPresence_biomes %>% 
  filter(!is.na(GROWTHFORM_STD)) %>% 
  group_by(biomes,GROWTHFORM_STD) %>% 
  summarise(N_species=n_distinct(Species))

biomes_sp<-spPresence_biomes %>% 
  filter(!is.na(GROWTHFORM_STD)) %>% 
  group_by(biomes) %>% 
  summarise(N_species=n_distinct(Species))

growth_biomes$biomes<-as.factor(growth_biomes$biomes)
growth_biomes$biomes<-factor(growth_biomes$biomes,levels = c("moist","tropical.mixed",
                                                             "savanna","grasslands",
                                                             "dry","xeric","mediterranean",
                                                             "temperate.mixed","coniferous","prairies","taiga","tundra"))

indx<-match(growth_biomes$biomes,biomes_sp$biomes)
growth_biomes$total_sp<-biomes_sp$N_species[indx]
growth_biomes$prop<-growth_biomes$N_species/growth_biomes$total_sp

growth_biomes_sub<-
  growth_biomes %>% 
  filter(!is.na(biomes))

pdf("./figs/Growth_composition_numbers.pdf")
ggplot(data = growth_biomes_sub, aes(x = "", y = prop, fill = GROWTHFORM_STD )) + 
geom_bar(stat = "identity", position = position_fill()) +
geom_text(aes(label = N_species), position = position_stack(vjust = 0.5), cex=2) +
coord_polar(theta = "y") +
facet_wrap(~ biomes)
dev.off()


pdf("./figs/Growth_composition_NO_numbers.pdf")
ggplot(data = growth_biomes_sub, aes(x = "", y = prop, fill = GROWTHFORM_STD )) + 
  geom_bar(stat = "identity", position = position_fill()) +
  coord_polar(theta = "y") +
  facet_wrap(~ biomes)
dev.off()

  