
library(phyloseq)
library(Maaslin2)
library(tidyverse)
library(ggsci)
library(reshape2)
library(ggsignif)
library(cowplot)
library(vegan)


setwd("/mnt/4tb/home/mgimenez/Matias/Agriculture/Agriculture_Figures")

physeq=readRDS("physeq_ITS.rds")

otu_tab <- otu_table(physeq)
cnm <- colnames(otu_tab)

ASVs <- character(0)
for (i in 1:length(cnm)){
  
  ASI <- paste0("ASV",i)
  ASVs <- c(ASVs, ASI)
}


taxa_names(physeq) <- ASVs

#Load metadata information
meta <- as.data.frame(sample_data(physeq))
tab <- read_delim("Metadata_Agricultura_16S.csv", delim = ",")
tab <- tab[-c(4,5),]

meta[,21] <- tab$Season
colnames(meta)[21] <- "Season"
meta[,22] <- as.numeric(tab$`Abundancia ind/m²`)
colnames(meta)[22] <- "EW_Abundance"
meta[,23] <- as.numeric(gsub(",", ".", tab$`Biomasa g/m²`))
colnames(meta)[23] <- "EW_biomass"
meta[,24] <- as.numeric(gsub(",", ".",tab$DAP))
colnames(meta)[24] <- "Aparent_density"
meta[,25] <- as.numeric(tab$`RH (%)`)
colnames(meta)[25] <- "RH"
colnames(meta)[2] <- "Land.use"
colnames(meta)[1] <- "Site"

sample_data(physeq) <- meta

LN1 <- gsub("CN", "NG", meta$Land.use)
meta$Land.use <- gsub("RA", "AR", LN1)


################################################################################
#Transform to even depth, filter only top 10 phyla
################################################################################

wh0 = genefilter_sample(physeq, filterfun_sample(function(x) x > 1), A=0.1*nsamples(physeq))
PG1 = prune_taxa(wh0, physeq)

#Transform even depth
GP1 = transform_sample_counts(PG1, function(x) 1E6 * x/sum(x))

#Keep only the 10 most abundant phyla
phylum.sum = tapply(taxa_sums(PG1), tax_table(PG1)[, "Phylum"], sum, na.rm=TRUE)
top10phyla = names(sort(phylum.sum, TRUE))[1:10]
PG1 = prune_taxa((tax_table(PG1)[, "Phylum"] %in% top10phyla), PG1)


spg <- as_tibble(sample_data(PG1))

 spg %>% 
  unite(treatment, c(Site, Land.use), sep = "_", remove = FALSE) -> spg1

 spg1$sample_names <- rownames(sample_data(PG1))
 spg2 <-  column_to_rownames(spg1, var="sample_names") 

 sample_data(PG1) <- spg2
 mergedGP = merge_samples(PG1, "treatment")


############################################
## Barplot Fungal community composition
############################################

physeq.phylum=tax_glom(mergedGP, taxrank="Phylum")

txt <- tax_table(mergedGP)
gns <- txt[,6]

otu_ms <- otu_table(physeq.phylum)
asv <- colnames(otu_ms)
txt <- tax_table(physeq.phylum)
rnt <- row.names(txt)
phy <- txt[,2]



idx <- which(rnt %in% asv)
phy.otu <- as.vector(phy[idx])

colnames(otu_ms) <- phy.otu

otu_phy_tab <- melt(otu_ms)
colnames(otu_phy_tab) <- c("Sample", "Phylum", "Rel_abundace")

phyn <- gsub("p__", "", otu_phy_tab$Phylum)

otu_phy_tab$Phylum <- phyn
otu_phy_tab %>% separate(Sample, c('Site', 'Land_use'), sep="_") -> otu_phy_tab

otu_phy_tab$Site <- gsub("Laurnaga", "L", otu_phy_tab$Site)
otu_phy_tab$Site <- gsub("Martirena", "M", otu_phy_tab$Site)
otu_phy_tab$Site <- gsub("Poni", "P", otu_phy_tab$Site)
otu_phy_tab$Site <- gsub("Salto", "SA", otu_phy_tab$Site)

otu_phy_tab$Land_use <- gsub("CN", "NG", otu_phy_tab$Land_use)
otu_phy_tab$Land_use <- gsub("RA", "AR", otu_phy_tab$Land_use)

# Stacked
ggplot(otu_phy_tab, aes(fill=Phylum, y=Rel_abundace, x=Land_use)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, size = 10),
        legend.key.height = unit(0.3, 'cm'),
        legend.position = c("top"),
        panel.background = element_blank())+
  facet_grid(~Site)+
  ylab("Relative Abundance")+
  xlab("Treatment")+
  ggtitle("A")+
  scale_fill_npg() -> comp_ITS

 comp_ITS

################################################################################
##Alpha diversity plot
################################################################################

df_rich <- as.data.frame(estimate_richness(PG1, split=TRUE, measures = "Observed"))

df_rich$Land.use <- as.factor(meta$Land.use)
df_rich$Site     <- meta$Site


alf_fun <- ggplot(data=df_rich, mapping = aes(x=as.factor(Land.use), y = Observed, color = Land.use))+
           geom_jitter(mapping=aes(shape=Site), size=2)+
           geom_boxplot(outlier.shape = NA, alpha = 0.1)+
           geom_signif(comparisons=list(c("NG", "AR")), annotations = "p.val<0.05")+
           theme_classic()+
           theme(legend.position = 'none')+
           xlab("Land use")+
           ylim(c(290, 610))+
           ggtitle(label = "B")+
           scale_color_aaas()
  
alf_fun

 ##############################
 ###Core microbiome analysis###
 ##############################


library(viridis)

PG1.rel <- microbiome::transform(PG1, "compositional")
core_16S <- microbiome::core(PG1.rel, detection=0.00001, prevalence=99/100)

cnm <- as.vector(tax_table(core_16S)[,6])
cnm1 <- make.unique(cnm)
cnm2 <- gsub("g__", "", cnm1)

taxa_names(core_16S) <- cnm2

prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)
detections <- format(round(detections, 4), nsmall = 4)
p1 <- microbiome::plot_core(core_16S, 
                plot.type = "heatmap", 
                colours = gray,
                prevalences = prevalences, 
                detections = detections, min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")

p1 <- p1 + theme_bw() + theme(axis.text.y=element_text(face = "italic")) + ylab("Fungal genus") + ggtitle("Fungal core taxa")
print(p1 + scale_fill_viridis())
tbl_ITS <- tax_table(core_ITS)

tiff("Supp_Fig2.tiff", units="in", width=7, height=4, res=600)
print(p1 + scale_fill_viridis())
dev.off()

################################################################################
#NMDS from phylum distribution
################################################################################

PG.ord <- ordinate(PG1, "NMDS", "bray")
p2 = plot_ordination(PG1, PG.ord, type="samples", color="Land.use", shape="Site") 
p2 + geom_polygon(aes(fill=Land.use), alpha = 0.5) + geom_point(size=3) + ggtitle("Samples")


meta1 <- meta[,c(4,6,12,13,18,23,24,25)]
cn1 <- colnames(meta1)
#cn1[4] <- "CEC"
cn1[3] <- "SOM"
cn1[5] <- "Clay"
colnames(meta1)<-cn1


env_nmds1 <- envfit(PG.ord, meta1, na.rm = TRUE)

plot(env_nmds1)

#Extract scores from nmds for ggplot2
data.scores = as.data.frame(scores(PG.ord)$site)
data.scores$Treatment <- meta$Land.use
data.scores$Site <- c(rep("L",4), rep("M",6), rep("P", 6), rep("SA",6))

#Multiply the scores
en_coord_cont = as.data.frame(scores(env_nmds1, "vectors")) # * ordiArrowMul(env_nmds)
en_coord_cat = as.data.frame(scores(env_nmds1, "factors")) #* ordiArrowMul(env_nmds)

rownames(en_coord_cont) <- c("pH", "Mg", "SOM", "P", "Clay", "EW_biomass", "AD", "RH%")

ggITS = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = Treatment, shape=Site), size = 3)+ 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, linewidth =0.5, colour = "grey50") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            label = row.names(en_coord_cont), size=3, hjust = -0.1, vjust=-0.1) + 
  stat_ellipse(aes(x=NMDS1,y=NMDS2,colour=Treatment),level = 0.9, alpha=0.5)+
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Treatment")+
  ggtitle(label = "C")+
  scale_color_aaas()

ggITS


# ASVs data frame to be explained 
tdf <-  as.data.frame(otu_table(PG1))
# Vector indicating land use for each sample
luc <- sample_data(PG1)$Land.use
# Vector explaining reptitions for sites and land use in order to constrain permutations
trt <- sample_data(PG1)$treatment

adon.results <-   adonis(tdf ~ luc,  method="bray",  perm=999)
#Adonis results: there is a significant difference in the separation of samples of diff. LUCs
adon.results$aov.tab

df <- adon.results$aov.tab$Df[1]
sum_of_squares <- adon.results$aov.tab$SumsOfSqs[1]
mean_squares <- adon.results$aov.tab$MeanSqs[1]
f_model <- adon.results$aov.tab$F.Model[1]
r2 <- adon.results$aov.tab$R2[1]
p_value <- adon.results$aov.tab$`Pr(>F)`[1]

plot_ord <- ggITS + annotate("text", x = Inf, y = Inf, label = paste("PERMANOVA p-value<", signif(p_value, 3)), hjust = 1.1, vjust = 1.1,color="gray35")
plot_ord1 <- plot_ord + annotate("text", x=Inf, y=Inf, label= "NMDS Stress= 0.104", hjust = 1.55, vjust = 2.3, color="gray35") 

ggITS <- plot_ord1
ggITS

###################################################
#### Maaslin2 differential abundance analysis  ####
###################################################


physeq.genus=tax_glom(physeq, taxrank="Genus")

otu_ms <- otu_table(physeq.genus)
asv <- colnames(otu_ms)

txt <- tax_table(physeq.genus)
rnt <- row.names(txt)
gns <- txt[,6]

idx <- which(rnt %in% asv)
gns.otu <- as.vector(gns[idx])

colnames(otu_ms) <- gns.otu

otu_df <- as.data.frame(otu_ms)
otu_df <- rownames_to_column(otu_df, "ID")
write_delim(otu_df, "OTU_table_ITS.tsv", delim="\t")
meta_df <- as.data.frame(meta)
meta_df1 <- cbind(rownames(meta_df),meta_df)
colnames(meta_df1)[1] <- "ID"
write_delim(meta_df1, "meta_table.tsv", delim="\t")


fit_data <- Maaslin2(
  "OTU_table_ITS.tsv", "meta_table.tsv", 'Output_ITS_pH_maaslin2',
  min_abundance = 1, min_prevalence = 0.1,
  fixed_effects = c('Land.use'), #P, Arena, Arcilla, Limo, CIC
  random_effects = c('pH', 'Site', 'Season'), #'CO'
  normalization = "CLR",
  plot_heatmap = TRUE,
  standardize = FALSE)

sig_tab <- read_delim("Output_ITS_pH_maaslin2/significant_results.tsv")

dif_tab <- sig_tab[which(sig_tab$metadata=="Land.use"),]


# change the factor levels so it will be displayed in correct order

ord <- order(dif_tab$coef)
dif_tab$feature <- factor(dif_tab$feature, levels = as.character(dif_tab$feature[ord]))

idx <- which(dif_tab$coef<0)
dif_tab$Class <- rep("Natural Grassland", dim(dif_tab)[1])
dif_tab$Class[idx] <- "Agriculture Rotation"

#df2 <- tab1[order(tab1$coef),]

dif_tab %>%
  mutate(coef1 = -(log(qval)*coef)) -> df2

dif_tab %>%
  mutate(log2FC = log2(exp(coef))) -> df2


df2$feature <- gsub("g__", "", df2$feature)
df2$feature <- gsub("g__", "", df2$feature)
df2$feature <- gsub("_gen_Incertae_sedis", "_gen_Inc_sed", df2$feature)

idx <- which(df2$qval<0.1)
df3 <- df2[idx,]


ord <- order(df3$log2FC)
df3$feature <- factor(df3$feature, levels = as.character(df3$feature[ord]))
colnames(df3)[10] <- "Land use"


ggplot(df3, aes(x = feature, y = coef)) +
  geom_bar(aes(fill = `Land use`), stat = 'identity') +  # color by class
  coord_flip() +  # horizontal bars
  geom_text(aes(y = 0, label = feature, hjust = as.numeric(log2FC > 0)), fontface="italic", size=3)+  # label text based on value
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none' ) +
  ylab("Effect size")+
  ylim(c(-15, 15))+
  ggtitle("D")+
  scale_fill_aaas()-> ef_size_plot

ef_size_plot

#plot_grid(ggITS, ef_size_plot, ncol = 2, rel_widths = c(2,1)) -> gr1
#plot_grid(comp_ITS, alf_fun, ncol = 2, rel_widths = c(2,1)) -> alf

plot_grid(alf_fun, ef_size_plot, nrow=2, rel_heights = c(1,3)) -> left_gr
          
plot_grid(comp_ITS, ggITS, nrow=2, rel_heights = c(1,1)) -> right_gr 
         

plot_grid(right_gr, left_gr, ncol=2, rel_widths = c(2,1)) -> Fig_1


Fig_1




tiff("Fig_1.tiff", units="in", width=13, height=8, res=800)
plot_grid(right_gr, left_gr, ncol=2, rel_widths = c(2,1)) -> Fig_1
Fig_1

dev.off()

###########################################
##Sampling sites and map
###########################################
#Supplementary figure


require(sf)
require(raster) 
require(maptools)

uy <- getData('GADM', country='Uruguay', level=1)


uy_sf <- st_as_sf(uy)


sitios<-read.csv("Sample_table.csv", header=TRUE,sep=",")

sitiosA=sitios[which(sitios$Treatment=="Agricultura"), ]
sitiosA$SampleID <- c("SA", "P", "L", "M")

sitios_sf <- st_as_sf(x = sitiosA, coords = c("Long", "Lat"))
st_crs(sitios_sf) <- "WGS84" # set the coordinate reference system

ggplot(data = uy_sf) +
    geom_sf(fill="darkseagreen1")+
  geom_sf(data = sitios_sf, mapping=aes(shape=SampleID), size=4)+
  theme(legend.title = element_blank()) -> map_uy


tiff("Supp_fig_1.tiff", units="in", width=10, height=7, res=800)

ggplot(data = uy_sf) +
  geom_sf(fill="darkseagreen1")+
  geom_sf(data = sitios_sf, mapping=aes(shape=SampleID), size=4)+
  theme(legend.title = element_blank()) -> map_uy
map_uy

dev.off()

########################################
### Analysis using fungal traits database (FUNGuild)
########################################
                              
tab <- read_delim("../matrix_count_FUNGuild.taxa.guilds.txt", delim = "\t")

Taxonomy <- txt[,6]
guild1 <- tab$trophicMode
guild <- gsub("Undephined", "Undefined", guild1)

# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(physeq.genus), "matrix")
# transpose if necessary
asv <- colnames(OTU1)

txt <- tax_table(physeq.genus)
rnt <- row.names(txt)
gns <- txt[,6]

idx <- which(rnt %in% asv)
gns.otu <- as.vector(gns[idx])
colnames(OTU1) <- gns.otu

OTU2 <- t(OTU1)
asv <- colnames(otu_ms)
OTUdf <- as.data.frame(OTU2)
OTUg <- cbind(OTUdf, guild1) 

#OTUtp <- cbind(OTUdf, tpmode) 

OTUg %>% 
  group_by(guild1) %>% 
  summarise(across(everything(), sum)) -> OTUguild

OTUg <- as.matrix(OTUguild[,2:23])
rownames(OTUg) <- OTUguild$guild1
meta <- sample_data(physeq.genus)
colnames(meta)[2] <- "Land use"

meta_df <- rownames_to_column(meta_df1)
colnames(meta_df)[1] <- "Samples"

OTUtp <- rownames_to_column(as.data.frame(OTUg))
colnames(OTUtp)[1] <- "TrophicMode"

write_delim(meta_df, "meta_table.tsv", delim="\t")
write_delim(OTUtp, "OTU_guild.tsv", delim="\t")

fit_data <- Maaslin2(
  "OTU_guild.tsv", "meta_table.tsv", 'Guild_maaslin2',
  min_abundance = 1, min_prevalence = 0.1,
  fixed_effects = c('Land.use'),
  random_effects = c('Site'),
  normalization = "CLR",
  plot_heatmap = TRUE,
  standardize = FALSE
)

tab <- read_delim('Guild_maaslin2/all_results.tsv', delim ="\t")

long_otu <- pivot_longer(OTUtp, cols= meta_df$Samples, names_to = "Samples")

unique(long_otu$TrophicMode)

idx <- which(long_otu$TrophicMode=="Pathotroph" | long_otu$TrophicMode=="Pathotroph-Saprotroph" | 
        long_otu$TrophicMode=="Saprotroph-Symbiotroph" |  long_otu$TrophicMode=="Symbiotroph")

otu_tp <- long_otu[idx,]
Site <- character(0)
treat <- character(0)

for(i in 1:length(otu_tp$Samples)){
  
  smp <- otu_tp$Samples[i]
  
  sit <- strsplit(smp, "-")[[1]][1]
  trt <- strsplit(smp, "-")[[1]][2]  
  
  Site <- c(Site, sit)
  treat <- c(treat, trt)
}

treat[which(treat=="CN1" | treat=="CN3" | treat=="CN2")] <- "NG"
treat[which(treat=="RA1" | treat=="RA3" | treat=="RA2")] <- "AR"

otu_tp$Site <- Site
otu_tp$Treatment <- treat

tpmode <- otu_tp$TrophicMode
tpm1 <- gsub("Pathotroph-Saprotroph", "Facultative pathotroph", tpmode)
tpm2 <- gsub("Pathotroph", "Obligate pathotroph", tpm1)
tpm3 <- gsub("Saprotroph-Symbiotroph", "Facultative symbiotroph", tpm2)
tpm4 <- gsub("Symbiotroph", "Obligate symbiotroph", tpm3)

otu_tp$TrophicMode <- tpm4

ggplot(data=otu_tp, mapping = aes(x=Treatment, y=value))+
  geom_violin(aes(fill=TrophicMode), alpha=.5)+
  geom_jitter(aes(shape=Site, color=TrophicMode,), width= 0.25, size=2.5)+
  stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
               colour = "black", width = 0.2)+
  theme_bw()+
  facet_wrap(~TrophicMode, scales = "free")+
  geom_signif(comparisons=list(c("NG", "AR")))+
  paletteer::scale_fill_paletteer_d("rcartocolor::Antique")+
  paletteer::scale_color_paletteer_d("rcartocolor::Antique") -> p4
p4 
             
tiff("Supp_Fig_4.tiff", units="in", width=13, height=8, res=800)
p4
dev.off()
