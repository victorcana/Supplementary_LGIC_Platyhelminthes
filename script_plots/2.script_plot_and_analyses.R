#setwd("/script_plots/")
getwd()
#install.packages("tidyverse")
library(tidyverse)
LGIC.todo <- read.delim("LGIC_Platyhelminthes.tsv",sep = "\t", header = TRUE) 
#LGIC.todo <- read.delim("LGIC_Platyheminthes_2.tsv",sep = "\t", header = TRUE) 


#################### Figure 2 ####################### #
##### Figure 2A: PCA script #####
#install_github("vqv/ggbiplot")
library(ggbiplot)

datos_nci <- LGIC.todo 
datos_nci <- datos_nci %>% mutate_at(vars(starts_with("Match")), ~ifelse(is.na(.), 1, .))%>% ###A los NA ponerles valor 1 para el análisis 
  filter(!str_detect(Class, "Free-living")) # Filter out free-living flatworms (No-Neodermata)
pca1 <- prcomp(datos_nci[, 2:7], scale = T)  # Perform PCA
summary(pca1) # Print PCA summary



ir.LGIC.family <- datos_nci[, 28] # LGIC family information
ir.class <- datos_nci[, 8] # Class taxonomic information
set.seed(123)
g <- ggbiplot::ggbiplot(pca1, obs.scale = 1, var.scale = 1,
                        groups = ir.LGIC.family, ellipse = TRUE, 
                        circle = FALSE, alpha = 0.0)+
  geom_point(aes(shape=ir.class, color=ir.LGIC.family), alpha = 0.3, size=2)
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')+ 
  scale_shape_discrete(name = "Neodermata class")+
  scale_colour_manual(values = c("#003f5c", "#58508d", "#009E73", "#ff6361","#ffa600"),name = "LGIC family")  + 
  scale_fill_manual(values = c("#003f5c", "#58508d", "#009E73", "#ff6361","#ffa600")) +
  scale_x_continuous(breaks = seq(-7, 5, 2), limits = (c(-7, 5))) +
  #scale_colour_viridis_d(alpha = .6,option = "C", direction = 1)+
  #  scale_fill_viridis_d( alpha = .6,option = "C", direction = 1)+
  theme_classic()

print(g)# Print PCA plot

# Save PCA plot (Figure 2A) in different formats
ggsave(plot = g, "plot_evalue_PCA.png", width = 22, height = 11,units = "cm",limitsize = FALSE)
ggsave(plot = g, "plot_evalue_PCA.svg",width = 22, height = 11,units = "cm", limitsize = FALSE)
ggsave(plot = g, "plot_evalue_PCA.pdf",width = 22, height = 11,units = "cm", limitsize = FALSE)





#### Figure 2B: K-means script #####
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/factoextra")
library(factoextra)
#install.packages("NbClust")
library("NbClust")
library(gridExtra)

# Modify data and filter out free-living flatworms
datos_nci <- LGIC.todo 
datos_nci <- datos_nci %>% mutate_at(vars(starts_with("Match")), ~ifelse(is.na(.), 1, .))%>% # Replace NA with 1 for analysis
  filter(!str_detect(Class, "Free-living"))  # Filter out free-living flatworms (No-Neodermata)

# Determine the number of clusters
set.seed(123)
nb <- NbClust(datos_nci[, 2:7], distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans")
## According to the majority rule, the best number of clusters is 2 and 4.
###The consequent analysis was performed by four clusters to obtain a higher partition number.

# Perform k-means clustering with 4 clusters
set.seed(123)
km.result <- kmeans(datos_nci[, 2:7], 4)  # Select 4 clusters
#km.result$cluster

# Save clustering results to a table (uncomment if needed)
#datos_ncikm <-cbind(datos_nci, Cluster = km.result$cluster)
#write.table(datos_ncikm, file = "datos_nci_clusteringkm.tsv", sep = "\t", row.names = FALSE)


# Evaluate each cluster through silhouette analysis
##Silhouette width is a metric that measures the quality of clustering.
#Provides an estimate of how well each data point clusters compared to how well it would fit in neighboring clusters.
##The silhouette width varies from -1 to 1, where higher values indicate well-defined and separated clusters
##An average silhouette width of 0.49 for a K-means clustering result can be considered a relatively good result,
km.sil <-cluster::silhouette(km.result$cluster, dist(datos_nci[, 2:7]))
set.seed(123)
fviz_silhouette(km.sil) +  
  scale_colour_manual(values = c("#003f5c", "#58508d", "#009E73", "#ff6361","#ffa600","#BCEE68"), name = "Cluster")  + 
  scale_fill_manual(values = c("#003f5c", "#58508d", "#009E73", "#ff6361","#ffa600","#BCEE68"),name = "Cluster") 

ggsave("plot_evalue_average_silhouete.png",width = 30, height = 15,units = "cm", limitsize = FALSE)
ggsave("plot_evalue_average_silhouete.pdf",width = 30, height = 15,units = "cm", limitsize = FALSE)
ggsave("plot_evalue_average_silhouete.svg",width = 30, height = 15,units = "cm", limitsize = FALSE)

# Plot Figure 2B with 4 selected clusters 
set.seed(123)
h <- fviz_cluster(km.result, data = datos_nci[, 2:7],labelsize = 0 , alpha=0.4,main =" ", obs.scale = 1, var.scale = 1,pointsize =2,
                  xlab="PC1 (82.7% explained var.)", ylab="PC1 (13.5% explained var.)")+
  geom_point( alpha = 0.0)+
  scale_shape_discrete(name = "Cluster                                        ")+
  scale_x_continuous(breaks = seq(-7, 5, 2), limits = (c(-7, 5))) +
  scale_colour_manual(values = c("#003f5c", "#58508d", "#009E73", "#ff6361","#ffa600","#BCEE68"), name = "Cluster                                        ")  + 
  scale_fill_manual(values = c("#003f5c", "#58508d", "#009E73", "#ff6361","#ffa600","#BCEE68"),name = "Cluster                                        ") +
  #  theme_ggstatsplot()
  theme_classic()
h
ggsave("plot_evalue_clustered_kmean.png",width = 22, height = 11,units = "cm", limitsize = FALSE)
ggsave("plot_evalue_clustered_kmean.pdf",width = 22, height = 11,units = "cm", limitsize = FALSE)
ggsave("plot_evalue_clustered_kmean.svg",width = 22, height = 11,units = "cm", limitsize = FALSE)



#### Plot combined Figure 2A and 2B #####
library(gridExtra)
aaa=list(g,h) # Combine PCA and K-means plots
plot.distribution <- grid.arrange(grobs = aaa,  ncol = 1) # Arrange plots in a single column

ggsave(plot  = plot.distribution,"plot_evalue_clustered_PCA_kmean.png",width = 22, height = 22,units = "cm", limitsize = FALSE)
ggsave(plot  = plot.distribution,"plot_evalue_clustered_PCA_kmean.svg",width = 22, height = 22,units = "cm", limitsize = FALSE)
ggsave(plot  = plot.distribution,"plot_evalue_clustered_PCA_kmean.pdf",width = 22, height = 22,units = "cm", limitsize = FALSE)


#### Figure 3: Distribution of E-values##### 
#For significant differences see Supplementary Figure S6
# install.packages("ggridges")
library(ggridges)
# install.packages("ggplot2")
library(ggplot2)
#install.packages("ggstatsplot")
library(ggstatsplot)
#if (!require(gridExtra)) {
#  install.packages("gridExtra")}
library(gridExtra)
#install.packages("svglite")
library(svglite) #Para guardar el gráfico con ggsave en formato svg

Evalue2_log10N_round_longer_Factor <- LGIC.todo %>%
  pivot_longer(cols = c(2,3,4,5,6,7), names_to = "Match", values_to = "Evalue") %>% 
  mutate(Match = factor(Match, levels = c("Match_Metazoa_evalue", "Match_Protostomia_evalue", "Match_Spiralia_evalue", "Match_Lophotrochozoa_evalue", "Match_Platyhelminthes_evalue", "Match_Neodermata_evalue")))

# Replace long names with short names
Evalue2_log10N_round_longer_Factor <- Evalue2_log10N_round_longer_Factor %>% 
  dplyr::mutate(Match = str_replace(Match, "Match_Metazoa_evalue", "Metazoa"),
         Match = str_replace(Match, "Match_Protostomia_evalue", "Protostomia"),
         Match = str_replace(Match, "Match_Spiralia_evalue", "Spiralia"),
         Match = str_replace(Match, "Match_Lophotrochozoa_evalue", "Lophotrochozoa"),
         Match = str_replace(Match, "Match_Platyhelminthes_evalue", "Platyhelminthes"),
         Match = str_replace(Match, "Match_Neodermata_evalue", "Neodermata")) %>%
  dplyr::mutate(Match = factor(Match, levels = c("Metazoa", "Protostomia", "Spiralia", "Lophotrochozoa", "Platyhelminthes", "Neodermata")))



# Generate plots for the distribution of E-values
Class_analysis <- Evalue2_log10N_round_longer_Factor %>% filter(!str_detect(Class, "Free-living")) %>% group_by(Class) %>% dplyr::summarise() %>% as.data.frame()
Family.LGIC <- Evalue2_log10N_round_longer_Factor %>%  group_by(LGIC.Family) %>% dplyr::summarise() %>% as.data.frame()
plot <- list()
for(j in 1:length(Family.LGIC[[1]])){
  for(i in 1:length(Class_analysis[[1]])){
    
    plot[[paste(Class_analysis[[1]][i],Family.LGIC[[1]][j])]] <- ggplot(Evalue2_log10N_round_longer_Factor %>% 
                                                                          filter(str_detect(Class, Class_analysis[[1]][i])) %>%
                                                                          filter(str_detect(LGIC.Family, Family.LGIC[[1]][j])), 
                                                                        aes(x = Evalue, y = Match, fill = stat(quantile) )) +
      stat_density_ridges(
        geom = "density_ridges_gradient", calc_ecdf = TRUE,
        quantiles = 4, quantile_lines = FALSE, scale = 1.5
      )+
      scale_fill_viridis_d(name = "Quartiles", alpha = .6,option = "D", direction = -1)+ 
      labs(x = "-log10(e-value)", y = "Taxa", subtitle = paste(Family.LGIC[[1]][j],"-", Class_analysis[[1]][i])) +
      scale_x_continuous(breaks = seq(0, 310, 50), limits = (c(0, 310))) +
      theme_ggstatsplot()
    #ggsave(paste("plot_",Family.LGIC,".png"))
  }
}

# Arrange plots into a grid
plot.distribution <- grid.arrange(grobs = plot,  ncol = 4)
grid.arrange(grobs = plot,  ncol = 4)

ggsave("plot_Distribution_evalue.png", plot=plot.distribution, width = 40, height = 20, units = "cm")
ggsave("plot_Distribution_evalue.pdf", plot=plot.distribution, width = 40, height = 20, units = "cm")
ggsave("plot.svg", plot=plot.distribution, width = 40, height = 20, units = "cm")



#### Supplementary Figure S6: Significant Differences ####
#install.packages("ggstatsplot")
library(ggstatsplot)

Evalue2_log10N_round_longer_Factor <- LGIC.todo %>%
  pivot_longer(cols = c(2,3,4,5,6,7), names_to = "Match", values_to = "Evalue") %>% 
  mutate(Match = factor(Match, levels = c("Match_Metazoa_evalue", "Match_Protostomia_evalue", "Match_Spiralia_evalue", "Match_Lophotrochozoa_evalue", "Match_Platyhelminthes_evalue", "Match_Neodermata_evalue")))

# Replace long names with short names
Evalue2_log10N_round_longer_Factor <- Evalue2_log10N_round_longer_Factor %>% 
  dplyr::mutate(Match = str_replace(Match, "Match_Metazoa_evalue", "Metazoa"),
         Match = str_replace(Match, "Match_Protostomia_evalue", "Protostomia"),
         Match = str_replace(Match, "Match_Spiralia_evalue", "Spiralia"),
         Match = str_replace(Match, "Match_Lophotrochozoa_evalue", "Lophotrochozoa"),
         Match = str_replace(Match, "Match_Platyhelminthes_evalue", "Platyhelminthes"),
         Match = str_replace(Match, "Match_Neodermata_evalue", "Neodermata")) %>%
  mutate(Match = factor(Match, levels = c("Metazoa", "Protostomia", "Spiralia", "Lophotrochozoa", "Platyhelminthes", "Neodermata")))

# Convert NA to 0 in E-values
Evalue2_log10N_round_longer_Factor1 <- Evalue2_log10N_round_longer_Factor %>%
  mutate(Evalue = ifelse(is.na(Evalue), 0, Evalue))

# Generate grouped comparisons of E-values
set.seed(123)
grouped_ggbetweenstats(
  data            = Evalue2_log10N_round_longer_Factor1 %>% 
    filter(!str_detect(Class, "Free-living")) %>%   
    arrange(LGIC.Family, Class,Species,Match,ID),
  #  ggplot.component = scale_y_continuous(limits = c(0, 310)),
  plotgrid.args = list(nrow = 4),
  x               = Match,
  y               = Evalue,
  type            = "np",
  xlab            = "Taxonomic level",
  ylab            = "-log10(e-value)",
  ggplot.component = list(theme(plot.title = element_text(size = 8, face = "bold")),
                          ###                          ggplot2::scale_y_continuous(sec.axis = ggplot2::dup_axis(name = "")), Quitar para el paper
                          theme(plot.subtitle = element_text(size = 8, face = "bold"))),
  grouping.var    = Class.LGIC,
)

# Save the plot in various formats
ggsave("plot_evalue_grouped_ggbetweenstats.png",width = 70, height = 50,units = "cm", limitsize = FALSE)
ggsave("plot_evalue_grouped_ggbetweenstats.pdf",width = 70, height = 50,units = "cm", limitsize = FALSE)
ggsave("plot_evalue_grouped_ggbetweenstats.svg",width = 70, height = 50,units = "cm", limitsize = FALSE)



#### Figure 4: Proportion of Proteins Without Orthologs ####
#install.packages("ggstatsplot")
library("ggstatsplot")

#Transform columns 2 to 7 into a single row.
Evalue2_log10N_round_longer_Factor <- LGIC.todo %>%
  pivot_longer(cols = c(2,3,4,5,6,7), names_to = "Match", values_to = "Evalue") %>% 
  mutate(Match = factor(Match, levels = c("Match_Metazoa_evalue", "Match_Protostomia_evalue", "Match_Spiralia_evalue", "Match_Lophotrochozoa_evalue", "Match_Platyhelminthes_evalue", "Match_Neodermata_evalue")))

# Replace long names with short names
Evalue2_log10N_round_longer_Factor <- Evalue2_log10N_round_longer_Factor %>% 
  dplyr::mutate(Match = str_replace(Match, "Match_Metazoa_evalue", "Metazoa"),
         Match = str_replace(Match, "Match_Protostomia_evalue", "Protostomia"),
         Match = str_replace(Match, "Match_Spiralia_evalue", "Spiralia"),
         Match = str_replace(Match, "Match_Lophotrochozoa_evalue", "Lophotrochozoa"),
         Match = str_replace(Match, "Match_Platyhelminthes_evalue", "Platyhelminthes"),
         Match = str_replace(Match, "Match_Neodermata_evalue", "Neodermata")) %>%
  mutate(Match = factor(Match, levels = c("Metazoa", "Protostomia", "Spiralia", "Lophotrochozoa", "Platyhelminthes", "Neodermata")))

# Calculate the number of proteins without orthologs
AA <- Evalue2_log10N_round_longer_Factor %>%
  dplyr::select(Match, Species, Evalue, Class, LGIC.Family) %>%
  filter(is.na(Evalue)) %>% group_by(Species,Match, LGIC.Family,Class) %>% dplyr::summarise(N_NA =n()) %>% ungroup()

# Calculate the total number of proteins 
BB <- Evalue2_log10N_round_longer_Factor %>%
  dplyr::select(Match, Species, Evalue, Class, LGIC.Family) %>%
  group_by(Species,Match, LGIC.Family,Class) %>% dplyr::summarise(N=n()) %>% ungroup()

# Combine the data
Numero_NA_Orden <-  full_join(AA,BB)%>%
  mutate(N_NA = replace(N_NA, is.na(N_NA), 0)) %>% 
  mutate(NA_prop = N_NA/N, Hit_prop = 1-NA_prop)%>%
  mutate(Especies_Match = paste(Match, Species, LGIC.Family, sep = "_")) %>% 
  mutate(Class.LGIC = paste(Class, LGIC.Family, sep = "_"))#%>% print(n=100)
#write.table(Numero_NA_Orden, file  = "Numero_NA_Orden", sep = "\t", row.names = FALSE)

# Assign identifier codes to each species
Especies  <- LGIC.todo  %>% group_by(Species)  %>% dplyr::summarise() 
Especies$subject  <- row.names(Especies)

# Join the data frames
Numero_NA_statis <- left_join(Numero_NA_Orden, Especies, by="Species")  
Numero_NA_statis <- Numero_NA_statis %>%
  select(last_col(), everything())


#### Normality Test ###
Datos_normalidad <- dplyr::filter(Numero_NA_statis,
                                  Class %in% c("Cestoda", "Trematoda", "Monopisthocotylea", "Polyopisthocotylea"),
                                  Match %in% c("Metazoa",  "Protostomia", "Spiralia","Lophotrochozoa", "Platyhelminthes","Neodermata")) %>%
  dplyr::select(Class, Match, NA_prop, Species,LGIC.Family)  %>% 
  pivot_wider(names_from = Match, values_from = NA_prop) 

shapiro.test(Datos_normalidad$Metazoa)
shapiro.test(Datos_normalidad$Protostomia)
shapiro.test(Datos_normalidad$Lophotrochozoa)
shapiro.test(Datos_normalidad$Platyhelminthes)
shapiro.test(Datos_normalidad$Neodermata)
shapiro.test(Datos_normalidad$Spiralia)
# Data is not normal


#### Significant Differences ###
# Generate grouped comparisons for Figure 4

set.seed(123)
Figure4 <- grouped_ggbetweenstats(
  data            = dplyr::filter(Numero_NA_statis,
                                  #                    Match %in% c("Lophotrochozoa",  "Neodermata"),
                                  Class %in% c("Cestoda", "Trematoda", "Monopisthocotylea", "Polyopisthocotylea")),
  
  x               = Match,
  y               = NA_prop,
  type            = "np",
  plotgrid.args = list(nrow = 2),
  ggplot.component = scale_y_continuous(limits = c(0, 1.5)),
  xlab            = "Taxonomic level", ##Pensar si quitar si se unen los gráficos
  ylab            = "Proportion of proteins without orthologs",
  grouping.var    = LGIC.Family
  #  pairwise.display = "all"
)
Figure4
# Save the plot in various formats
ggsave(plot = Figure4, "Figure4_plot_NA_ggbetweenstats.png", width = 45, height = 30,units = "cm",limitsize = FALSE)
ggsave(plot = Figure4, "Figure4_plot_NA_ggbetweenstats.pdf",width = 45, height = 30,units = "cm", limitsize = FALSE)
ggsave(plot = Figure4, "Figure4_plot_NA_ggbetweenstats.svg",width = 45, height = 30,units = "cm", limitsize = FALSE)