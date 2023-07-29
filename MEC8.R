# https://rpubs.com/angela34a

# Package loading ####
library(FEAST) # for source tracking
library(tidyverse)
library(vegan)
library(FSA) # for dunn post-hoc test
library(caret) 
library(broom)
library(patchwork) # for combining plots into 1
library(ggfortify)
library(pairwiseAdonis) # for permanova post-hoc




# Load data ####

# load metadata
metadata <- read.csv("~/faks/Ecology and Ecosystems/MEC-8 (Specific research project)/data/metadata_Lunz_revised.csv", 
                     row.names = 1)

# load asv table
asv_table <- read.csv("~/faks/Ecology and Ecosystems/MEC-8 (Specific research project)/data/asv_table.csv", row.names=1)

# identifier need to be identical between asv_table and metadata
identical(rownames(metadata), colnames(asv_table))
# if there is some that are different use setdiff(rownames(metadata), colnames(asv_table)) to see which

# data on sources/sinks
feast_metadata <- read.csv("~/faks/Ecology and Ecosystems/MEC-8 (Specific research project)/data/feast_metadata.csv", row.names=1)
# identifiers for sources have to be unique - source 1 needs to be Source1_1, Source1_2, ... )
# SampleID column need to be  rownames


# FEAST ####
# classes of datasets

# feast_metadata should be a data.frame
class(feast_metadata)

#asv should be matrix and each column an ineger in order to perform feast
feast_asv <- as.matrix(asv_table)                # matrix as a whole
for(i in 1:ncol(feast_asv)){                 # each column as integer
  feast_asv[, i] <- as.integer(feast_asv[, i])
}


# whole  table should be a matrix
class(feast_asv)

# individual column should be integer
class(feast_asv[, 1])



# final step - form needs to be as with vegan (sites rows, species columns)
feast_asv <- t(feast_asv)
feast_asv[1:10, 1:5]


#to remove from environment
rm(i) 


# perform FEAST
feast_output <- FEAST(C = feast_asv,              # transposed asv_table
                      metadata = feast_metadata,  # sink or source information
                      different_sources_flag = 0, # 0 if sources are not assigned to different sinks
                      EM_iterations = 1000,       # number of iterations
                      dir_path = "C:/Users/Angela Cukusic/Documents/faks/Ecology and Ecosystems/MEC-8 (Specific research project)/Results", 
                      outfile = "feast_results")  # prefix for saving output file


#we got an output, import it now
feast_table <- read.table("C:/Users/Angela Cukusic/Documents/faks/Ecology and Ecosystems/MEC-8 (Specific research project)/Results/feast_results_source_contributions_matrix.txt", 
                          sep = "\t", row.names = 1)


## Raw contribution to plottable data ####
# feast_output is of class "list"
# we want to make it a matrix

# 1. make a new vector with source types (and it is "character")
sources <- paste("source", c(1:8), sep = "")


# 2. make empty matrix with nrow = number of SINKS (5), 
# and ncol = number of possible sources (20+ samples are designated to 1-8 sources)
source_contributions <- matrix(NA,                        # empty
                               nrow = nrow(feast_table),  # 5 sinks
                               ncol = length(sources))    # n sources

# 3. fill matrix 
for(i in 1:length(sources)){    # for each source add contribution 
  source_contributions[, i] <- rowSums(feast_table[,      # sum up contributions for all 5 rows
                                                   grepl(sources[i], colnames(feast_table))]) 
  # of only the columns which have sourceX in name
}

# 4. we added contributions assigned to each of 8 sources, now add from unknown sources
source_contributions <- cbind(source_contributions, feast_table$Unknown)


# 5. check if all contributuons for individual sink amount to 1
rowSums(source_contributions) 


# 6. convert to data frame and bring back rownames and colnames
source_contributions <- as.data.frame(source_contributions)
colnames(source_contributions) <- c(sources, "unknown") # col
source_contributions <- source_contributions %>%        # row
  mutate(sample_id = rownames(feast_table), .before = 1) 


#  final clean up of environment
rm( sources, i)


## Plot source contribution ####


### all sinks ####
source_contributions %>% 
  
  pivot_longer(cols = -sample_id, names_to = "sources", values_to = "contributions") %>% 
  mutate(sources = as.factor(sources))     %>% 
  
  
  # find  relative abundances for each date in each water kind (well/pump/surf)
  group_by(sample_id) %>%
  mutate(rel_cont = contributions / sum(contributions) * 100) %>%
  ungroup() %>%
  
  # plotting
  ggplot(aes(x = sample_id, y = rel_cont, fill = sources)) +
  geom_bar(stat = "identity") +
  
  # labelling   
  geom_text(aes(label =  
                  # label if % is big enough
                  ifelse(rel_cont > 2, 
                         # if yes then  paste, if no then empty ""
                         paste0(round(rel_cont, 4) , "%"), "" )  ),
            # how to position it
            position = position_stack(vjust = 0.5), size = 4) +
  
  # theme  
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  theme_classic() + 
  scale_fill_brewer(palette = "Set1") +
  labs(x= "Sinks", y="Relative contributions %", fill = "Sources") +
  scale_x_discrete(labels = c("Untersee \n(31m)", 
                              "Untersee \n(5m)", 
                              "Untersee - \nnear OSB inflow", 
                              "Untersee - \nat OSB \ninflow", 
                              "Sediment: \nUntersee - \nat OSB \ninflow")) 



### inflow sinks ####
feast_table %>% 
  rownames_to_column("sample_id") %>% 
  # find the two sinks in question  
  filter(sample_id %in% c("OSB_in_LUS_inflow_LUS_at_OSB_inflow",
                          "OSB_in_LUS_sed_LUS_sediment")) %>% 
  
  # find the individual samples from source 6
  select( "sample_id",  contains("source6")) %>% 
  # remove the last part of names, which says "source6"  
  rename_all(~str_remove(., ".{10}$")) %>% 
  
  pivot_longer(cols = -sample_id, names_to = "sources", values_to = "contributions") %>% 
  mutate(sources = as.factor(sources))     %>% 
  
  
  # find  relative abundances for each date in each water kind (well/pump/surf)
  group_by(sample_id) %>%
  mutate(rel_cont = contributions / sum(contributions) * 100) %>%
  ungroup() %>%
  
  # plotting
  
  ggplot(aes(x = sample_id, y = rel_cont, fill = sources)) +
  geom_bar(stat = "identity", color = "black") +
  
  # labelling   
  geom_text(aes(label =  
                  # label if % is big enough
                  ifelse(rel_cont > 2, 
                         # if yes then  paste, if no then empty ""
                         paste0(round(rel_cont, 4) , "%"), "" )  ),
            # how to position it
            position = position_stack(vjust = 0.5), size = 4) +
  
  # theme  
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  theme_classic() + 
  labs(x= "Sinks", y="Relative contributions %", fill = "Sources") +
  scale_x_discrete(labels = c("Untersee - \nat OSB \ninflow", 
                              "Sediment: \nUntersee - \nat OSB \ninflow")) +
  scale_fill_manual(values = c( "#006400", "#FFD700", "lightblue", "#FFD700",
                                "lightblue", "#006400", "lightblue" ))


# Sample type distribution ####

metadata_type <- metadata %>% 
  rownames_to_column("sources") %>% 
  select("sources", "type_III") %>% 
  mutate(type_III = case_when(
    # assign categories  
    type_III %in% c("Creek", "Hyporheic",
                    "Pond", "Lake", "Inflow") ~ "Surface", 
    TRUE ~ type_III)) 



feast_table %>% 
  rownames_to_column("sink_id") %>% 
  # rename the columns so that they can be connected to the metadata (type info)  
  rename_all(~str_remove(., ".{10}$")) %>% 
  pivot_longer(cols = - sink_id, names_to = "sources", values_to = "contributions") %>% 
  # filter out contributions from unknow sources
  filter(sources != "Unknown") %>% 
  # join with metadata
  left_join(., metadata_type, by = "sources") %>% 
  filter(type_III != "NA") %>% 
  filter( contributions < 0.05) %>% 
  # plot
  ggplot(aes(x= type_III, y= contributions)) + 
  geom_boxplot(aes(fill = sink_id) ) +
  theme_linedraw() +
  scale_fill_manual(labels=c("Untersee (31m)", 
                             "Untersee (5m)",
                             "Untersee (near inflow)", 
                             "Untersee (at inflow)", 
                             "Sediment (at inflow)"), 
                    values=c("#acc7d9", "#8bb0ca","#5180a2","#3e647d","#d8b365")) +
  theme(legend.position=c(0.3,0.8), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14), 
        legend.background = element_rect(), 
        legend.box.background  =   element_rect(color = "gray")) +
  labs(y = "Contributions", x= "Sample type") + 
  scale_x_discrete(labels = c("Biofilm", "Ground water", "Sediment", 
                              "Soil", "Spring water", "Surface water") )



# Spatial distribution ####

metadata_rank <- metadata %>% 
  rownames_to_column("sources") %>% 
  select("sources", "rank") %>% 
  filter(rank != "sink")

feast_table %>% 
  rownames_to_column("sink_id") %>% 
  # rename the columns so that they can be connected to the metadata (rank info)  
  rename_all(~str_remove(., ".{10}$")) %>% 
  pivot_longer(cols = - sink_id, names_to = "sources", values_to = "contributions") %>% 
  # filter out contributions from unknow sources
  filter(sources != "Unknown") %>% 
  # join with metadata
  left_join(., metadata_rank, by = "sources") %>% 
  filter(rank != "NA") %>% 
  # plot
  ggplot(aes(x= rank, y= contributions)) + 
  geom_boxplot( fill = "violet") +
  facet_wrap(~sink_id, scales = "free_y",  labeller = 
               as_labeller( c("LUS_dp_31m_LUS_31m" = "Untersee (31m)",
                              "LUS_dp_5m_LUS_5m" = "Untersee (5m)",
                              "LUS_in_OSB_LUS_near_OSB_inflow" = "Untersee (near inflow)",
                              "OSB_in_LUS_inflow_LUS_at_OSB_inflow" = "Untersee (at inflow)",
                              "OSB_in_LUS_sed_LUS_sediment" = "Sediment (at inflow)" ))) +
  theme_linedraw() + 
  labs(y = "Contributions", x= "Ranked by spatial distance ")




# Pre-processing for community analysis ####

# number of reads per sample
asv_table %>% colSums() %>% summary() 

# data has  been rarefied by Lucas before forwwarding





# Alpha diversity ####

shannondiv <- vegan::diversity(t(asv_table))
# recalculating it as the effective number of species to be more comparable across studies.
ens <- exp(shannondiv) %>% as.data.frame() %>% rownames_to_column("sample_id") %>% 
  rename("diversity" = ".")

diversity_data <- specnumber(t(asv_table)) %>% as.data.frame() %>% rownames_to_column("sample_id") %>% 
  rename("richness" = ".") %>% inner_join(., ens, by ="sample_id" )

# we will plot diversity by source so make clean categorization here
feast_metadata %>% 
  # remove last two strings
  mutate(Env = str_remove(Env, ".{2}$")) %>% 
  # where SourceSink column is a source keep the value from the Env column
  mutate(Source_category = case_when(SourceSink == "Source" ~ Env,
                                     # otherwise keep the value from SourceSink column                                   
                                     TRUE ~ SourceSink) ) %>% 
  # some trouble with source 3
  mutate(Source_category = case_when(Source_category %in% c("source3", "source3_") ~ "source3",
                                     TRUE ~Source_category)) %>% 
  
  rownames_to_column("sample_id")  %>% 
  select("sample_id", "Source_category") %>% 
  mutate(Source_category = as.factor(Source_category)) %>% 
  # we now have clean category and can connect it to the diversity data
  inner_join(., diversity_data,   by = "sample_id") %>% 
  select(-sample_id) %>% 
  
  #pivot longer so we can facet
  pivot_longer(cols = -Source_category, names_to = "measure", values_to = "values") %>% 
  
  # plot
  ggplot( aes(x = Source_category, y = values)) +
  #stat_boxplot(geom ='errorbar') + 
  geom_boxplot(fill = "violet") +
  theme_linedraw()+
  coord_flip() + 
  facet_wrap(~measure, scales = "free",  labeller = 
               as_labeller( c("diversity" = "Shannon's diversity index", 
                              "richness" = "Richness") ) )+
  theme(legend.position = "none",
        axis.title=element_text(size=14,face="bold"),
        axis.text.y = element_text(size=10,face="bold"),
        strip.text = element_text(size = 14)) +
  labs( x = "Source category",
        y = NULL) +
  scale_x_discrete(labels = c("Sink",  "Source 1", "Source 2", "Source 3", "Source 4", 
                              "Source 5", "Source 6",  "Source 7",  "Source 8"))  



## Test it ####

# 1. the same data manipulation as above
feast_metadata %>% 
  mutate(Env = str_remove(Env, ".{2}$")) %>% 
  mutate(Source_category = case_when(SourceSink == "Source" ~ Env,                           
                                     TRUE ~ SourceSink) ) %>% 
  mutate(Source_category = case_when(Source_category %in% c("source3", "source3_") ~ "source3",
                                     TRUE ~Source_category)) %>% 
  
  rownames_to_column("sample_id")  %>% 
  select("sample_id", "Source_category") %>% 
  mutate(Source_category = as.factor(Source_category)) %>% 
  inner_join(., diversity_data,   by = "sample_id") %>% 
  select(-sample_id) -> data_for_alpha


# 2. perform ANOVA/ kruskal-wallis

# 1. assumption of ANOVA: Normal distribution: 
# Each population to be compared should be normally distributed 
# (NOT THE ENTIRE DATASET).
shapiro.test (  data_for_alpha[which(data_for_alpha$Source_category == "source1"), "diversity"] ) # NOT ND
shapiro.test (  data_for_alpha[which(data_for_alpha$Source_category == "source2"), "diversity"] ) # ND
shapiro.test (  data_for_alpha[which(data_for_alpha$Source_category == "source3"), "diversity"] ) # ND
shapiro.test (  data_for_alpha[which(data_for_alpha$Source_category == "source4"), "diversity"] ) # NOT ND
shapiro.test (  data_for_alpha[which(data_for_alpha$Source_category == "source5"), "diversity"] ) # NOT ND
shapiro.test (  data_for_alpha[which(data_for_alpha$Source_category == "source6"), "diversity"] ) # ND
# source 7 is too small sample size to perform the test  
shapiro.test (  data_for_alpha[which(data_for_alpha$Source_category == "source8"), "diversity"] ) # ND
shapiro.test (  data_for_alpha[which(data_for_alpha$Source_category == "Sink"), "diversity"] ) # ND

# A: would be better to do a non-parametric test


## 2. assumption of ANOVA: Homogeneity of variance: 
## Variance in the populations compared should be the same/similar. 
bartlett.test ( diversity ~ Source_category, data= data_for_alpha )
# assumption met


# 3. test
kruskal.test (diversity ~ Source_category, data= data_for_alpha) 
# 4. post-hoc
FSA::dunnTest(diversity ~ Source_category, data= data_for_alpha)



# Beta diversity ####
bc_diss <- vegdist(t(asv_table), method = "bray")
nmds <- metaMDS(bc_diss)


# make palette for plotting
pal = c("green4", "navy", "orange", "brown", "blue", "lightblue")


# lets plot an nmds where the color is sample type, shape is source sink
as.data.frame(nmds$points) %>% 
  cbind(data_for_alpha$Source_category, metadata_type$type_III) %>% 
  rename("source_cat" = "data_for_alpha$Source_category",
         "type" = "metadata_type$type_III") %>% 
  mutate(source_cat = case_when(source_cat == "Sink" ~ "sink",
                                TRUE ~ "source")) ->  data_for_nmds
  
# plot 
data_for_nmds %>%   
  ggplot(aes(x = MDS1, y = MDS2)) +
  geom_point(aes(fill = type, color = type, 
                 shape= source_cat), 
             size = 5, stroke=1, alpha = 0.6) +
# customizing the plot  
  scale_fill_manual(values =  pal) + 
  scale_color_manual(values = pal) +
  scale_shape_manual(values = c(21,23)) +
  
  theme_linedraw() +
  annotate(geom="text",x=-1.4, y=1.5, label="stress =", color="black") +   
  annotate(geom="text", x=-1.1, y=1.5, label=round(nmds$stress,4),
           color="black") +
  
  #stat_ellipse(aes(colour = type), linewidth = 0.3) + 
  #stat_ellipse(aes(fill = type), geom="polygon",level=0.95,alpha=0.04) +
  
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "gray")) +
  labs(color = "Water \ntype", shape = "Source") +
  guides(fill = FALSE) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) -> plot_nmds
  

plot_nmds




## Test community differences ####

# test differences in group centroid with adonis (anova version for multivariate statistics)
permanova <- adonis2(t(asv_table) ~ data_for_nmds$type, 
                    method = "bray" , permutations = 999)
permanova # 0.001
pair_res <- pairwiseAdonis::pairwise.adonis2(t(asv_table) ~ type, data = data_for_nmds)
# Spring_vs_Surface, GW_vs_Spring are only one not significantly distinct!!


# test dispersion as an assumption (variance test for multivariate statitics)
disp <- betadisper( vegdist(  t(asv_table), method = "bray" ), data_for_nmds$type) 
anova(disp)
TukeyHSD(disp) 
# Surface-GW have diff dispersions, does not interfere with centroid differences





# Environmental data ####

## pca (environmenta variables of water communities) ####

df_pca <- metadata %>% 
  # we only have measurements for water, not for soil nor sediment 
  filter(type_I == "Water") %>% 
  select(where(is.numeric)) %>% 
  select(-c("filter", # dont need it
            "N_NO2", "P_PO4")) # NAs
 

# remove the highly correlated ones, like internalATP and totalATP
cor_matrix <- cor(df_pca)
highly_correlated <- caret::findCorrelation(cor_matrix, cutoff = 0.9)
df_pca <- df_pca[, -highly_correlated]

# now there is 22 observations and 10 environmental variables
# lets see if any env vars can be packed neatly in PCs
pca_env <- prcomp(df_pca, scale = TRUE)
biplot(pca_env)

# i believe it is beneficial to remove the sample PX_gw since it is a single outlier in almost every variable, and stresses the ordination
pca_env <- prcomp(df_pca[rownames(df_pca) != "PX_gw",], scale = TRUE)


# we need the column for coloring by type, get it from the original metadata
 metadata[rownames(df_pca) != "PX_gw", c(colnames(df_pca), "type_III", "type_I")] %>% 
   # remove the outlier sample, the columns that are highly correlated and soil/sed samples
  filter(type_I == "Water") -> data_for_pca

autoplot(pca_env, data =data_for_pca , 
         colour = 'type_III',
         size=5, loadings=TRUE, 
         loadings.colour="black", 
         loadings.label=TRUE, loadings.label.size=5, 
         loadings.label.colour="black" )+ 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(shape= 21, color = "black", size = 5) + 
  theme_linedraw() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  scale_color_brewer(palette = "Accent")+
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "gray"))  -> pca_plot

pca_plot


# Extract the loadings from the PCA result and convert to tidy format 
loadings_df <-   tidy(pca_env, matrix = "rotation")


# PC1
# Create the ggplot2 object for plotting the loadings
loadings_df %>% filter(PC == "1") %>% 
  #x=reorder(class,-amount,sum)
ggplot( aes(x = reorder(column, -value), y = value, group = factor(column))) +
  geom_bar(stat = "identity", position = "dodge", 
           fill= "gray90", color = "gray20") +
  geom_text(aes(label = column), position = position_dodge(width = 0.9), 
            #vjust = "inward",  
            angle = 90 , hjust = "inward") +
  labs(x = "Variables of PC1", y = "Loadings") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "none")  -> plot_pca1

# PC2 
loadings_df %>% filter(PC == "2") %>% 
  #x=reorder(class,-amount,sum)
ggplot( aes(x = reorder(column, -value), y = value, group = factor(column))) +
  geom_bar(stat = "identity", position = "dodge", 
           fill= "gray90", color = "gray20") +
  geom_text(aes(label = column), position = position_dodge(width = 0.9), 
            #vjust = "inward",  
            angle = 90 , hjust = "inward") +
  labs(x = "Variables of PC2", y = "Loadings") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "none") -> plot_pca2

# PC3
loadings_df %>% filter(PC == "3") %>% 
  #x=reorder(class,-amount,sum)
ggplot( aes(x = reorder(column, -value), y = value, group = factor(column))) +
  geom_bar(stat = "identity", position = "dodge", 
           fill= "gray90", color = "gray20") +
  geom_text(aes(label = column), position = position_dodge(width = 0.9), 
            #vjust = "inward",  
            angle = 90 , hjust = "inward") +
  labs(x = "Variables of PC3", y = "Loadings") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "none") -> plot_pca3

plot_pca1 + plot_pca2 + plot_pca3



# Test pca dispersions

# we can check the significance in difference between groups in both pca and nmds

euclidean_dist <- data_for_pca %>% select(-c("type_I" , "type_III")) %>% 
                  scale() %>% as.data.frame()  %>% 
                  dist(., method = "euclidean") 

disp_pca <- betadisper( euclidean_dist, data_for_pca$type_III )
anova(disp_pca) #no sign. differences 


disp_pca$distances %>% as.data.frame() %>% 
  rename("Distance" = ".") %>% 
  rownames_to_column("sample_id") %>% 
  cbind(data_for_pca) %>% 
# plot
  ggplot( aes(x = type_III, y = Distance)) +
  geom_boxplot( fill = "violet") +
  labs(x = NULL, y = "Environmental differences \nstandardized Euclidean distance") +
  theme_linedraw() +
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) -> pca_disp


## nmds (water communities) ####

asv_water <- asv_table [,rownames(df_pca)] %>% select(!"PX_gw")
bc_diss_water <- vegdist(t(asv_water), method = "bray")
nmds_water <- metaMDS(bc_diss_water)


as.data.frame(nmds_water$points) %>% 
  cbind(data_for_pca$type_III) %>% 
  rename("type" = "data_for_pca$type_III") %>% 
  
# plot 
  ggplot(aes(x = MDS1, y = MDS2)) +
  geom_point(aes(fill = type), shape = 21, 
             size = 5, color= "black" ) +
# customizing the plot  
  scale_fill_brewer(palette = "Accent") +
  theme_linedraw() +
  annotate(geom="text",x=-1.4, y=1.5, label="stress =", color="black") +   
  annotate(geom="text", x=-1.1, y=1.5, label=round(nmds_water$stress,4),
           color="black") +
  
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "gray")) +
  labs(color = "Water \ntype", shape = "Source") +
  guides(fill = FALSE) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) -> nmds_plot_water
  


euclidean_dist <- data_for_pca %>% select(-c("type_I" , "type_III")) %>% 
                  scale() %>% as.data.frame()  %>% 
                  dist(., method = "euclidean") 

disp_nmds <- betadisper( bc_diss_water, data_for_pca$type_III )
anova(disp_nmds) # sign. differences 


disp_nmds$distances %>% as.data.frame() %>% 
  rename("Distance" = ".") %>% 
  rownames_to_column("sample_id") %>% 
  cbind(data_for_pca) %>% 
# plot
  ggplot( aes(x = type_III, y = Distance)) +
  geom_boxplot( fill = "violet") +
  labs(x = NULL, y = "Bray-curtis distance") +
  theme_linedraw() +
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) -> nmds_water_disp


 nmds_plot_water  + pca_plot +
  nmds_water_disp + pca_disp 



## Constrained ordiantion ####


data_for_pca %>% select(where(is.numeric)) %>% 
  scale() %>% as.data.frame() -> env.data

decostand(t(asv_water), method = "hellinger") -> sp.dat

# check
identical(rownames(env.data), colnames(asv_water)) #TRUE

# perform model
simpleRDA <- capscale(sp.dat ~  ., data=env.data ,
                      distance = "bray")

vif.cca(simpleRDA) # all variables lower than 10


# test ####

# test env paameters
anova.model3 <- anova(simpleRDA, step=1000, by = "term")
anova.model3
# sign: bix* , fi***, TCC***, Cond, O2**,  Temp*

# new model 
simpleRDA <- capscale(sp.dat ~  bix + fi + TCC + Cond + O2 + Temp, data=env.data ,
                      distance = "bray")
# test env paameters
anova.model3 <- anova(simpleRDA, step=1000, by = "term")

# test of all canonical axes
anova.model <- anova.cca(simpleRDA, by='axis', step=1000)  # CAP1,CAP2,CAP3 significant 

# test the whole model
anova.model2 <- anova.cca(simpleRDA, step=1000) # sign
RsquareAdj(simpleRDA)$adj.r.squared # explains 30%



# plot ####

# vectors
ccavectors <- as.matrix(scores(simpleRDA, display = "bp", scaling = "sites")*3.32) %>% 
  t() %>% as.data.frame() %>% 
  rename("bix*" = "bix",
         "Temp*" = "Temp",
         "O2**" = "O2",
         "fi***" = "fi",
         "TCC***" = "TCC") %>%
  t() %>% as.data.frame()

#bix* , fi***, TCC***, Cond, O2**,  Temp*


# site coordinates
site_data <- scores(simpleRDA, display = "sites") %>% 
  as.data.frame() %>% 
  cbind(., data_for_pca$type_III) %>% 
  rename("type" = "data_for_pca$type_III")

# plotting
plot_cca <- 
  site_data %>% 
  ggplot( aes(x = CAP1, y = CAP2)) +
  geom_point(aes(fill = type), shape = 21, 
             size = 5, color= "black" )+  
  geom_segment(data = ccavectors, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               size = 1.2, arrow = arrow(length = unit(0.5, "cm"))) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_text(data = ccavectors, aes(x = CAP1*1.2, y = CAP2*1.1, 
                                   label = rownames(ccavectors)),
            size=6 ) +
  theme_linedraw() +
  scale_fill_brewer(palette = "Accent") +
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "gray")) 

plot_cca


