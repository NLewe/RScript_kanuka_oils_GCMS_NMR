# R Script for 
# "Comparison of chemical profiles of Kānuka (Kunzea robusta, Myrtaceae) essential oils"
# authors: Natascha Lewe,1,† Michaela Young,1,† Jan Vorster,2 Bella Paenga,3 Damian Skinner,3 Nikki Harcourt,3,4 Peter de Lange,5 Tia Haira,1 Storm Blockley-Powell,1 Andrew Munkacsi,1,* Robert Keyzers2,*
  
# Please cite the packages that are used.
# Package citation can be found in R with 
# citation ("tidyverse") etc.


# Install and load packages ####
# All packages can be installed with:
# install.packages ("tidyverse")

# Load package in R ##

library (tidyverse)  # basic package for data wrangling
library (readxl)     # to import excel data
library (FactoMineR) # for easy visualisation of PCA results 
library (factoextra) # calculation and visualisation of PCA
library (vegan)      # statistical test
library(ggpubr)     # helper package for ggplot - visualisation
library (ggrepel) # helper package for ggplot - visualisation



#  Load data into R ####
# All data and meta data is in excel tables.

## NMR data ####
# Preparation of the data:
# The NMR data was binned as described in the paper and the resulting values copied into an excel sheet.
# 
# 

NMR <-  read_excel ("data/NMR_data.xlsx", sheet = "K_oil_peak") 
# the NMR data is read into R from an excel file that is 
#located in the folder "data" in the working directory of R.
# Only the sheet "K_oil_peak" is loaded into R. 
# Colomns in the excel sheet are named b bin, rows are named with the sample ID

meta_NMR  <- read_excel ("data/NMR_data.xlsx", sheet = "meta_NMR")
# an excel sheet that has the sample ID  as rowname and consist of colums of meta data, for example origin of the sample (named grup here)

## NMR data clean up ####
# To remove noise and keep valuable bins, all bins with area values smaller than 0.000025 were removed #

NMR_test <- 
  meta_NMR %>%  
  select (sampleID, Sample, group) %>% 
  left_join(NMR)  %>%                    # unites meta data and NMR data
  select_if (~any (.>= 0.000025)) %>%   # removes small values
  select (!sampleID) %>%  # clean up of data table
  select (!Class) %>% 
  data.frame(row.names = "Sample") # change to a data frame to run PCA

## First tests of PCA results####
PCA_NMR <- PCA (NMR_test [,1:178], scale.unit = T, quali.sup = 1, graph = T, ncp = 5) 


## Check PCA calculation####

# Vidualise a scree plot to check how many dimensions of the PCA are relevant #

eig.val <- PCA_NMR$eig  # get the eigenvalues from PCA results
barplot(eig.val[, 2], # plot ar barplot (scree plot)
        names.arg = 1:nrow(eig.val),
        main = "Variances Explained by PCs (%)",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.val), eig.val[, 2],
      type = "b", pch = 19, col = "red")

# Usually, the number of dimension (principal componenets) needed to explain about 80 % of the variance in the data is a good cut-off

# PCA Plot for NMR ####
# Prepare lots in ggplot ##

# get data points for plot #
PCA_metric_data   <- PCA_NMR$ind$coord  %>%  as_tibble(rownames = "Sample") %>% 
  left_join(meta_NMR %>%  select (group, Sample) %>%  unique ())

# get variance explained by PC 1PC2 and PC3
PCA_eig_m_Dim1 <- round (PCA_NMR$eig[1,2],1)
PCA_eig_m_Dim2 <- round (PCA_NMR$eig[2,2],1)
PCA_eig_m_Dim3 <- round (PCA_NMR$eig[3,2],1)
# total variance explained by first three PCs #
PCA_eig_m_Dim1 + PCA_eig_m_Dim2 +PCA_eig_m_Dim3

#define colors
cols =  c("Commercial samples"=   "#999999", "Horomaka"= "#009E73", "Taranaki" = "#0072B2", 
          "Tairāwhiti"= "#D55E00", "Taitokerau"= "#CC79A7")



# Plot NMR PCA ###
plotPca  <- 
  PCA_metric_data %>%  
  rename ("Origin" = "group") %>% 
  ggplot (aes (x = Dim.1, y = Dim.2, color = Origin,label = Sample)) + 
  geom_point (size =5) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_text_repel(aes (), show.legend = F, point.size = 5) +
  theme_bw() + 
  xlab(label = "PC1 (32.4 %)") +
  ylab ("PC2 (17.8 %)") + 
  theme (legend.position = "bottom", legend.text = element_text(size = 14), legend.title = element_text(size = 15),
         axis.text = element_text(size = 11), axis.title = element_text (size = 12)) +
  scale_color_manual (values =cols)

# dimension 2, 3
plotPca2  <-
  PCA_metric_data %>%  
  rename ("Origin" = "group") %>% 
  ggplot (aes (x = Dim.2, y = Dim.3, color = Origin, label = Sample)) + 
  geom_point (size =5) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_text_repel(aes (), show.legend = F, point.size = 5) +
  theme_bw() + 
  xlab(label = "PC1 (17.8 %)") +
  ylab ("PC2 (11.5 %)") + 
  theme ( legend.text = element_text(size = 14), legend.title = element_text(size = 15),
         axis.text = element_text(size = 11), axis.title = element_text (size = 12)) +
  scale_color_manual (values = cols)

#ggarrange (plotPca2,plotPca2, nrow = 2, labels = c("A", "B") )


# PCA NMR contributions ####
NMR_cont1 <- fviz_contrib(PCA_NMR, choice = "var", axes = 1 , top = 50) +
  theme_bw() +
  ggtitle("Contributions of NMR bins to PC1") +
  xlab ("bin") +
  theme (axis.text.x = element_text(angle =90, size = 8))

NMR_cont2 <- fviz_contrib(PCA_NMR, choice = "var", axes = 2, top = 20) +
  theme_bw() +
  ggtitle("Contributions of NMR bins to PC2") +  
  xlab ("bin") +
  theme (axis.text.x = element_text(angle =90, size =8))



NMR_cont3 <- fviz_contrib(PCA_NMR, choice = "var", axes = 3, top = 20) +
  theme_bw() +
  ggtitle("Contributions of NMR bins to PC3") + 
  xlab ("bin") + 
  theme (axis.text.x = element_text(angle =90, size =8))


ggarrange (NMR_cont1, ggarrange(NMR_cont2, NMR_cont3, nrow = 1, ncol = 2,labels = c("B", "C"), widths =c(1,1) ),
           labels = c("A", NULL), nrow = 2)






#### GC ####

GC_data <-  read_excel ("data/Kanuka_RAK.xlsx", sheet = "Kanuka_list")

meta_GC  <- read_excel ("data/kanuka_RAK.xlsx", sheet = "meta_GC")

### tidy GC data, are the NAs =  zeroes?

GCMS_data <- 
  GC_data [,8:20] %>% 
  as_tibble () %>%  
  select (!Sample_M1) %>% 
  pivot_longer(!Name, names_to =  "sampleID" ) %>%  
  mutate (value = replace_na (value,0)) %>%  
  pivot_wider(names_from = Name, values_from = value) %>%
  left_join (meta_GC) %>% 
  data.frame(row.names = "sampleID", check.names = F) 



# calculate PCA ####
PCA_GCMS <- PCA (GCMS_data [,1:74], scale.unit = T, quali.sup = 74, graph = T, ncp = 5)




PCA_metric_data_GC   <- PCA_GCMS$ind$coord  %>%  
  as_tibble(rownames = "sampleID") %>% 
  left_join(meta_GC %>%  select (group, Sample, sampleID) %>%  unique ())

PCA_eig_m_Dim1_GC <- round (PCA_GCMS$eig[1,2],1)

PCA_eig_m_Dim2_GC <- round (PCA_GCMS$eig[2,2],1)
PCA_eig_m_Dim3_GC <- round (PCA_GCMS$eig[3,2],1)
PCA_eig_m_Dim1_GC + PCA_eig_m_Dim2_GC +  PCA_eig_m_Dim3_GC



#Plot GC PCA ######
plotPca_GC  <- 
  PCA_metric_data_GC %>%  
  rename ("Origin" = "group") %>% 
  ggplot (aes (x = Dim.1, y = Dim.2, color = Origin, label = Sample)) + 
  geom_point (size =5) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_text_repel(aes (), show.legend = F, point.size = 5) +
  theme_bw() + 
  xlab(label = "PC1 (27.6 %)") +
  ylab ("PC2 (18.0 %)") + 
  theme (legend.position = "bottom", legend.text = element_text(size = 14), legend.title = element_text(size = 15),
         axis.text = element_text(size = 11), axis.title = element_text (size = 12)) +
  scale_color_manual (values = cols)

plotPca_GC2  <- 
  PCA_metric_data_GC %>%  
  rename ("Origin" = "group") %>% 
  ggplot (aes (x = Dim.2, y = Dim.3, color = Origin, label = Sample)) + 
  geom_point (size =5) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_text_repel(aes (), show.legend = F, point.size = 5) +
  theme_bw() + 
  xlab(label = "PC2 (18 %)") +
  ylab ("PC3 (12.6 %)") + 
  theme ( legend.text = element_text(size = 14), legend.title = element_text(size = 15),
         axis.text = element_text(size = 11), axis.title = element_text (size = 12)) +
  scale_color_manual (values = cols)

#colorblind palette
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## PC 1 for NMR and GCMS  PLOT 
ggarrange (plotPca, plotPca_GC, common.legend = T, ncol = 2, legend = "bottom",  
           labels = c("A) NMR", "B) GCMS" ))

# scree plot to check how many dimensions of the PCA are relevant # 

eig.val <- PCA_GCMS$eig
barplot(eig.val[, 2], 
        names.arg = 1:nrow(eig.val), 
        main = "Variances Explained by PCs (%)",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.val), eig.val[, 2], 
      type = "b", pch = 19, col = "red")



## Plot PCA tests # 

plot.PCA (PCA_GCMS, axes = c(1,2), habillage = 74, title = "PCA of GC/MS results")  # habillage is the group infromation, the gravitational center is shown with the group name

plot.PCA (PCA_GCMS, axes = c(2,3), habillage = 74, title = "PCA of GC/MS results")

GC_cont1 <- 
  fviz_contrib (PCA_GCMS, choice = "var", axes = 1, top = 50) +
  theme_bw() +
  ggtitle("Contributions to PC1") +
  xlab (NULL) +
  theme (axis.text.x = element_text(angle =45, hjust = 1, size = 8))

GC_cont2 <- 
  fviz_contrib (PCA_GCMS, choice = "var", axes = 2, top = 20) +
  theme_bw() +
  ggtitle("Contributions to PC2") +
  xlab (NULL) +
  theme (axis.text.x = element_text(angle =45, hjust = 1, size = 8))

GC_cont3 <- 
  fviz_contrib (PCA_GCMS, choice = "var", axes = 3, top = 20) +
  theme_bw() +
  ggtitle("Contributions to PC3") +
  xlab (NULL) +
  theme (axis.text.x = element_text(angle =45, hjust = 1, size = 8))


## gg
ggarrange (plotPca2,NMR_cont1, nrow = 2, labels = c("A", "B") )

ggarrange (GC_cont1, ggarrange(GC_cont2, GC_cont3,align = "h", nrow = 1, ncol = 2,labels = c("B", "C"), widths =c(1,1) ),
           labels = c("A", NULL), nrow = 2)

# Plot PCA GC only 7 compound s of interest ####

PCA_GC_7 <- GCNM_cor  %>%  select (compound) %>%  unique () %>%  left_join( (GC_data [,8:20] %>%    as_tibble () %>%   select (!Sample_M1)), by = c("compound"= "Name") ) %>%  pivot_longer(!compound, names_to =  "sampleID" ) %>%  
  mutate (value = replace_na (value,0)) %>%  
  pivot_wider(names_from = compound, values_from = value) %>%
  left_join (meta_GC) %>% 
  data.frame(row.names = "sampleID") 



PCA_GCMS_7 <- PCA (PCA_GC_7 [,1:8], scale.unit = T, quali.sup = 8, graph = T, ncp = 3)




PCA_metric_data_GC_7   <- PCA_GCMS_7$ind$coord  %>%  
  as_tibble(rownames = "sampleID") %>% 
  left_join(meta_GC %>%  select (group, Sample, sampleID) %>%  unique ())

PCA_eig_m_Dim1_GC7 <- round (PCA_GCMS_7$eig[1,2],1)

PCA_eig_m_Dim2_GC7 <- round (PCA_GCMS_7$eig[2,2],1)
PCA_eig_m_Dim3_GC7 <- round (PCA_GCMS_7$eig[3,2],1)
PCA_eig_m_Dim1_GC7+ PCA_eig_m_Dim2_GC7 +  PCA_eig_m_Dim3_GC7

PCA_arrows_metrics<-
  PCA_GCMS_7$var$coord %>%  
  as_tibble (rownames = "metric") %>% 
  dplyr::rename("D1end" = "Dim.1", "D2end"= "Dim.2")



#PLOT ##
plotPca_GC_7  <- 
  PCA_metric_data_GC_7 %>%  
  rename ("Origin" = "group") %>% 
  ggplot (aes (x = Dim.1, y = Dim.2, color = Origin, label = Sample)) + 
  geom_point (size =5) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_text_repel(aes (), show.legend = F, point.size = 5) +  
  geom_segment(data = PCA_arrows_metrics, aes (x=0, xend= D1end*1, y = 0, yend = D2end*1), 
               arrow = arrow(length = unit(0.3, "picas")), color = "darkblue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics, aes (x  = D1end*1, y = D2end*1, label = metric), 
                    color = "darkblue", inherit.aes = F , force = 0.6) + 
  theme_bw() + 
  xlab(label = "PC1 (48.2 %)") +
  ylab ("PC2 (21.3 %)") + 
  theme (legend.position = "bottom", legend.text = element_text(size = 14), legend.title = element_text(size = 15),
         axis.text = element_text(size = 11), axis.title = element_text (size = 12)) +
  scale_color_manual (values = c( "#999999", "#009E73",  "#0072B2", "#D55E00", "#CC79A7"))

GC_7_cont1<- fviz_contrib(PCA_GCMS_7, choice = "var", axes = 1) +
  theme_bw () +
  ggtitle("Contributions to PC1")

GC_7_cont2<- fviz_contrib(PCA_GCMS_7, choice = "var", axes = 2) +
  theme_bw () +
  ggtitle("Contributions to PC2")

GC_7_cont3<- fviz_contrib(PCA_GCMS_7, choice = "var", axes = 3) +
  theme_bw () +
  ggtitle("Contributions to PC3")


ggarrange (plotPca_GC_7, (ggarrange (GC_7_cont1, GC_7_cont2, nrow = 1, ncol = 2, labels = c("B", "C"))),
           nrow = 2, ncol = 1, heights = c(1.7, 1), labels = c("A", NULL))

plotPca_GC2_7  <- 
  PCA_metric_data_GC_7 %>%  
  rename ("Origin" = "group") %>% 
  ggplot (aes (x = Dim.2, y = Dim.3, color = Origin, label = Sample)) + 
  geom_point (size =5) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_text_repel(aes (), show.legend = F, point.size = 5) +
  theme_bw() + 
  xlab(label = "PC2 (21.3 %)") +
  ylab ("PC3 (18.3 %)") + 
  theme (legend.position = "bottom", legend.text = element_text(size = 14), legend.title = element_text(size = 15),
         axis.text = element_text(size = 11), axis.title = element_text (size = 12)) +
  scale_color_manual (values = c( "#999999", "#009E73",  "#0072B2", "#D55E00", "#CC79A7"))




## PC 1 for NMR and GCMS  PLOT 
ggarrange (plotPca, plotPca_GC, common.legend = T, ncol = 2, legend = "bottom",  
           labels = c("A) NMR", "B) GCMS" ))

# scree plot to check how many dimensions of the PCA are relevant # 

eig.val <- PCA_GCMS$eig
barplot(eig.val[, 2], 
        names.arg = 1:nrow(eig.val), 
        main = "Variances Explained by PCs (%)",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.val), eig.val[, 2], 
      type = "b", pch = 19, col = "red")



## Plot PCA tests # 

plot.PCA (PCA_GCMS, axes = c(1,2), habillage = 74, title = "PCA of GC/MS results")  # habillage is the group infromation, the gravitational center is shown with the group name

plot.PCA (PCA_GCMS, axes = c(2,3), habillage = 74, title = "PCA of GC/MS results")

GC_cont1 <- fviz_contrib (PCA_GCMS, choice = "var", axes = 1, top = 20, title = "Contribution of variable to PC1")
GC_cont2 <- fviz_contrib (PCA_GCMS, choice = "var", axes = 2, top = 10, title = "Contribution of variable to PC2")
GC_cont3 <- fviz_contrib (PCA_GCMS, choice = "var", axes = 3, top = 10, title = "Contribution of variable to PC3")




# correlation coefficients####

GCNM_cor <- read_excel("data/Kanuka_RAK.xlsx", sheet = "corrCof")

corrCoef  <- 
  GCNM_cor %>%  
  pivot_wider( values_from = "value_perc", names_from = "NG") %>%  
  group_by (compound) %>% 
  summarize (cor = cor (NMR, GCMS, method = "pearson"))


corrPvalue <- 
  GCNM_cor %>%  
  pivot_wider( values_from = "value_perc", names_from = "NG") %>%  
  group_by (compound) %>% 
  summarize (corP = cor.test (NMR, GCMS, method = "pearson", alternative = "two.sided")$p.value)
##Procrustes analysis between PCAs#####


## Needs 
# symmetry = T - both solutions are scaled to#
# unit variance
# results in procrustes m2



procr <- procrustes (PCA_NMR$ind$coord, PCA_GCMS$ind$coord, symmetric = T)
#choice = c(1,2)  )  # with choice the number of dimensions is chosen
summary (procr)


# Plotting procrustes results##

plot (procr, kind =2
)
#Kind 2 plots show the residuals for each sample. 
#This allows identification of samples with the worst fit. 
#The horizontal lines, from bottom to top, are 
#the 25% (dashed), 50% (solid), and 75% (dashed) quantiles of the residuals.


##Test of Significance

#Function protest is a permutational test of the significance of the procrustes result.
#It is based on the correlation from a symmetric Procrustes analysis.


procrtest <- protest(X = PCA_NMR$ind$coord , Y = PCA_GCMS$ind$coord , 
                     scores = "sites", permutations = 99999)








#### NMDS test #### 

NMDS_NMR_center   <- metaMDS (NMR_center_test [,2:79], distance = "bray", trymax=1000, k =5, trace=F) #trace = F , then it stays quiet, k is number of dimension
#NMDS more robust than methods that include the magnitude of the distances, 
# get scores
scrs <- scores (NMDS_NMR_center, display = "sites")
scrs <-as.data.frame(scrs) %>%  as_tibble(rownames = "sampleID") %>%  left_join (meta_NMR)
center <- aggregate(cbind (NMDS1, NMDS2) ~KaMa1, data = scrs, FUN= mean) %>% 
  dplyr::rename("cNMDS1" = "NMDS1", "cNMDS2" = "NMDS2") %>% 
  left_join (meta_NMR %>% unique ()) 


segs <- left_join (scrs, center) 


#Plot NMDS in ggplot using phyloseq
# col  <- brewer.pal(5,  "Accent") # use all 8 colors from  palette
# new_pal  <- colorRampPalette(col)  #then use with col = new_pal (21) for 21 different colors


p <- 
  ggplot (segs , aes (x = NMDS1, y = NMDS2, color = KaMa1))  + 
  geom_segment(segs, mapping = aes (xend = cNMDS1, yend = cNMDS2)) + 
  geom_point(data = segs, mapping =aes (x = NMDS1, y= NMDS2), size = 3) + 
  geom_point (data = center , aes (x= cNMDS1, y = cNMDS2), size = 5)



p  + 
  guides (color = guide_legend (order =2,  title = "Plant species", 
                                label.theme = element_text(face = "italic", family = "sans", size = 11),ncol = 1), 
          shape = guide_legend(order=1, title = "Plant family",legend.position = "bottom" ))  +
  #scale_color_manual(values = new_pal(21))  + 
  theme_bw() + 
  # stat_ellipse(level=0.95) +
  scale_shape_manual (values= c(16,17, 18, 3,12,8 )) +
  annotate ("text", x= 0.2, y=0.4, label = "Stress = 0.11") #annotation is added by hand here, see NMDS_dm$stress
# NMDS_dm$stress
#the stress value reflects how well the ordination summarizes the observed distances among the samples.
#stress value >0.2 poor and might not be interpretable, values <0.1 good, between 0.1 and 0.2 usable but see dexter et all - stress depends on samples etc, 
#some of the distances might be misleading

##Test ecological null model
# result<-oecosimu(comm = df, method = "quasiswap", 
#                  nestfun = metaMDS, autotransform = FALSE, k = 5, #add dimension
#                  distance = "raup", nsimul = 1000,  statistic = "stress", #add distance measure
#                  alternative = "less", trace = FALSE, maxit = 200, 
#                  trymax = 50, sratmax = 0.9999999)

# result$oecosimu$statistic # stress of my data NMDS
# result$oecosimu$means
# result$oecosimu$z
# result$oecosimu$pval #this is relevant for reporting

# Similar to a permutation test  - in oecosimu() the values of the statistic are computed using data generated from one of the Null models
# How MANY of the stress values generated using the Null model of the data are less than or equal to the observed stress value (OR: alternative = "two.sided" (=different from), alternative = "greater"(or greater than or equal to))   .
# Plot the stress values of simulated NMDS
#Plot quick and dirty to check
# hist(as.vector(result$oecosimu$simulated), xlim = c(0,max(result$oecosimu$simulated)+.05), xlab = "Stress value",ylab = "Frequency observed", main = "", breaks = 7)
# abline(v = result$oecosimu$statistic, col = "red", lty = 2) 


# #Pretty ggplot  stress values randomised vs stress observed
# simVector<-as.data.frame (as.vector(result$oecosimu$simulated))
# simVector<-cbind(simVector,seq(1:1000))
# 
# stresshisto <- ggplot(simVector,aes(x=as.vector(result$oecosimu$simulated))) +  
#   geom_histogram(col="black",size=.25) +
#   theme_minimal() + 
#   xlim(c(0.07, .2)) + ylim(c(0,600)) + #change limits as fitting to stress value
#   geom_segment(aes(x = result$oecosimu$statistic, y = 0, xend = result$oecosimu$statistic), yend = 600,linetype=2,color = "red",size=.5) + #change limits of plot as needed
#   labs(x="Stress value")  +
#   labs(y="Frequency observed")
# 

pair.perm_raup_core <- read_rds ("resultsCh2/pair_perm_raup_core.rds")
permdisp_pairs_raup_core <- read_rds ("resultsCh2/permdisp_pairs_raup_core.rds")
permanova_raup_core <- read_rds ("resultsCh2/permanova_raup_core.rds")
permdisp_raup_core <- read_rds ("resultsCh2/permdisp_raup_core.rds")





