
library(tidyverse)
library(Maaslin2)
library(edgeR)
library(phyloseq)
library(ggsci)
library(ggsignif)
library(ggrepel)
library(paletteer)

#############################
## Nitrogen cycle Analysis ##
#############################

N_tab <- read_delim("FPKM_Ncyc_genes.tsv", delim = "\t")

meta <- read_delim("Metadata_Agricultura_16S.csv", delim = ",")
idx <-  grep(1, meta$'Muestra ')

meta_n <- read_delim("meta_new.tsv", delim = ",")
meta <- meta_n[idx,]
#meta$Sample  <-  smps

colnames(meta)[1] <- "Sitio"
colnames(meta)[2] <- "Treatment"
meta1 <- as.data.frame(meta)
rownames(meta1) <- colnames(N_tab)[-1]

N_FPKM <- column_to_rownames(N_tab, var = "target")

fit_data <- Maaslin2(
  N_FPKM, meta1, 'Out_NCyc2_maaslin2',
  min_abundance = 0.005,
  fixed_effects = c('Treatment'), #P, Arena, Arcilla, Limo, CIC
  #random_effects = c('Sitio'), #'CO'
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  # standardize = FALSE,
  normalization = 'CLR'
)

ntb <- read_delim(file="Out_NCyc2_maaslin2/all_results.tsv")

ntb %>%
  mutate(FC=exp(coef),
         log2FC=log2(FC),
         log10FC=log10(FC),
         log10qval=-log10(qval),
         log10pval=-log10(pval)) -> ntb1 

ntb1$Land_use <- rep("NO", dim(ntb1)[1])
ntb1$Land_use[which(ntb1$qval<0.25 & ntb1$coef<0)] <-"NG"
ntb1$Land_use[which(ntb1$qval<0.25 & ntb1$coef>0)] <-"AR"


ntb1$dalabel <- NA
ntb1$dalabel[ntb1$Diff_Abundance != "NO"] <- ntb1$feature[ntb1$Diff_Abundance != "NO"]


 ntb2<-ntb1[ntb1$qval<0.14,]
 
 ord <- order(ntb2$log2FC)
 ntb2$feature <- factor(ntb2$feature, levels = as.character(ntb2$feature[ord]))
 
 ntb2$Land_use <- rep("NG", length(ntb2$feature))
 ntb2$Land_use[ntb2$coef>0] <- "RA" 
 
 prc <- c("Denitrification","Assimilatory nitrate reduction","Denitrification", "Org. N Metabolism",
   "Org. N Metabolism","Dissimilatory nitrate reduction","Denitrification",rep("Org. N Metabolism",3),
   "Nitrogen fixation","Dissimilatory nitrate reduction", "Org. N Metabolism", "Assimilatory nitrate reduction",
   "Assimilatory nitrate reduction", "Denitrification", "Org. N Metabolism", "Org. N Metabolism",
   "Dissimilatory nitrate reduction","Org. N Metabolism")
 
 ntb2$Process <- prc
 

###########################
### Ncyc genes barplot ####
###########################

   
   lvln <- c("napA","nosZ","nirK","nirS", "nasA","nasB","nrfC","nirD","nirB","ureC", "gdh_K00262","glsA","gs_K00266","gs_K00265",
             "gdh_K15371","nmo","ureB","nifH","nirA",
             "gdh_K00261")
   ntb2$feature <- factor(ntb2$feature, levels = as.character(lvln))  
    
   ggplot(ntb2, aes(x = feature, y = coef)) +
     geom_bar(aes(fill = prc), alpha=.8, stat = 'identity') +  # color by class
     coord_flip() +  # horizontal bars
     geom_text(aes(y = 0, label = feature, hjust = as.numeric(log10FC > 0)), fontface="italic", size=4.5)+  # label text based on value
     geom_text(aes(x= 1, y=-.5, label="NG"), fontface="plain", family="sans", size=4, color= "gray50")+
     geom_text(aes(x= 1, y=.17, label="AR"), fontface="plain", family="sans", size=4, color= "gray50")+
      theme(axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),
           axis.title.y = element_blank(),
           panel.background = element_blank(),
           legend.title = element_blank(),
           legend.key.size = unit(0.2, 'cm'), #change legend key size
           legend.key.height = unit(0.2, 'cm'), #change legend key height
           legend.key.width = unit(0.2, 'cm'), #change legend key width
           legend.text = element_text(size=9),
           legend.position = c(.3,.8))+
     ylim(c(-0.5, 0.17))+
     ylab("Effect size")+
     ggtitle("D")+
     scale_fill_paletteer_d("MoMAColors::Clay")-> ncyc_plot
   ncyc_plot  
   
#########################
### Pcycle_metabolism ###
#########################
   
   
   
   #Read genes table
   PCtab <- read_delim('Pcyc_genes_database.tsv')
   #Read FPKM table
   FPKM <- read_delim("Pcyc_RotFPKM.tsv", delim="\t")
   FPKM <- FPKM[,c(1:26)]
   
   #Convert metagenome genes to enzyme names
   
   vles <- numeric(0)
   
   for(j in 1:length(FPKM$target)){
     
     gni <- FPKM$target[j]
     match <- grep(gni, PCtab$Gene_catalogue, fixed = TRUE)
     vles <- c(vles, match)
   }
   
   FPKM$target <- PCtab$Enzyme[vles]
   
   #Group FPKM counts by genes
   FPKM %>% 
     group_by(target) %>%
     summarise(across(everything(), sum),
               .groups = 'drop') %>%
     as.data.frame() -> FPKtab_group
   
   #Select only Agriculture samples and reshape df
   smps <- c("L-CN.Prd", "L-RA.Prd", "M-CN.Prd", "M-RA.Prd",  "P-CN.Prd",
             "P-RA.Prd", "SA-RA.Prd", "SA-CN.Prd")  
   
   idx <- which(colnames(FPKtab_group) %in% smps)
   FPK_df <- FPKtab_group[,c(1,idx)]
   cnms <- gsub(".Prd", "", colnames(FPK_df)) 
   colnames(FPK_df) <- cnms
   FPK_df <-  FPK_df[1:117,]
   rownames(FPK_df) <- FPKtab_group$target[1:117]
   
   P_FPKM <- FPK_df[,2:9]
   
   fit_data <- Maaslin2(
     P_FPKM, meta1, 'Out_PCyc2_maaslin2',
     min_abundance = 0.005,
     fixed_effects = c('Treatment'), #P, Arena, Arcilla, Limo, CIC
     #random_effects = c('Sitio'), #'CO'
     plot_heatmap = TRUE,
     plot_scatter = TRUE,
     # standardize = FALSE,
     normalization = 'CLR'
   )
   
   ptb <- read_delim(file="Out_PCyc2_maaslin2/all_results.tsv")
   
   ptb %>%
     mutate(FC=exp(coef),
            log2FC=log2(FC),
            log10FC=log10(FC),
            log10qval=-log10(qval),
            log10pval=-log10(pval)) -> ptb1 
   
   ptb1$Land_use <- rep("NO", dim(ptb1)[1])
   ptb1$Land_use[which(ptb1$qval<0.25 & ptb1$coef<0)] <-"NG"
   ptb1$Land_use[which(ptb1$qval<0.25 & ptb1$coef>0)] <-"AR"
   
   ptb1$dalabel <- NA
   ptb1$dalabel[ptb1$Diff_Abundance != "NO"] <- ptb1$feature[ptb1$Diff_Abundance != "NO"]
   
   ptb2<-ptb1[ptb1$qval<0.14,]
   
   ord <- order(ptb2$log2FC)
   ptb2$feature <- factor(ptb2$feature, levels = as.character(ptb2$feature[ord]))
   
   ptb2$Diff_Abundance <- rep("NG", length(ptb2$feature))
   ptb2$Diff_Abundance[ptb2$coef>0] <- "AR" 
   
   
   
   ##########################
   ### Pcyc genes barplot ###
   ##########################

   library(readxl)
   
   Proc_tab <- read_excel("40168_2022_1292_MOESM1_ESM.xlsx")
   Pr_tab <- Proc_tab[2:142,]
   Process <- character()
   
   for(i in 1:length(ptb2$feature)){
     
     tbg <- ptb2$feature[i]
     idx <- which(Pr_tab$Gene%in%tbg)[1]
     Proc <- Pr_tab$`Metabolic processes`[idx]
     Process <- c(Process, Proc)
     
   }
   
   Process1 <- gsub("and phosphinate", "",Process)
   Process2 <- gsub("Pyrimidine metabolism", "Org. P metabolism",Process1)
   Process3 <- gsub("Purine metabolism", "Org. P metabolism",Process2)
   ptb2$Process <- gsub("Pyruvate metabolism", "Org. P metabolism",Process3)
   
   lvl2 <- c("purT", "purB","ppdK","purH","ppc","purS","dut","purL", "pyrE","purC",
             "pstS","phnT","ugpC","pstB","phnC", "prsA",
             "gdh", "ppk","phnK","phoB",
             "phoH","pstC", "tmk","nrdD","dcd","pps","ndk")
   
   ptb2$feature <- factor(ptb2$feature, 
                          levels = as.character(lvl2))
   
   ggplot(ptb2, aes(x = feature, y = coef, group=Process))+
     geom_bar(aes(fill = Process), stat = 'identity', alpha=.9)+
     coord_flip() +  # horizontal bars
     geom_text(aes(y = 0, label = feature, hjust = as.numeric(log10FC > 0)), fontface="italic", size=3.5)+  # label text based on value
     geom_text(aes(x= 1, y=-1, label="NG"), fontface="plain", family="sans", size=4, color= "gray50")+
     geom_text(aes(x= 1, y=1, label="AR"), fontface="plain", family="sans", size=4, color= "gray50")+
     theme(axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),
           axis.title.y = element_blank(),
           panel.background = element_blank(),
           legend.title = element_blank(),
           legend.key.size = unit(0.2, 'cm'), #change legend key size
           legend.key.height = unit(0.2, 'cm'), #change legend key height
           legend.key.width = unit(0.2, 'cm'), #change legend key width
           legend.text = element_text(size=9),
           legend.position = c(.7,.6))+
     ylim(c(-1, 1.6))+
     ylab("Effect size")+
     ggtitle("E")+
     scale_fill_paletteer_d("rcartocolor::Antique") -> pcyc_plot
  
     pcyc_plot  
       
   
  
  ###########################################
  ### Cazymes Differential abundance
  ###########################################
   
   FPKtab <- read_delim("RotFPKM.tsv", delim = "\t")
   tax <- read_delim("CAZy_taxonomy.tsv", delim = "\t")
   
   ndx <- numeric(0)
   
   for(i in 1:length(FPKtab$target)){
     
     tgi <- FPKtab$target[i]
     idx <- which(tax$Gene_catalogue == tgi)
     ndx <- c(ndx, idx)
     
   }
   
   caz <- tax$CAZyFamily[ndx]
   
   FPKtab$target <- caz
   
   FPKtab %>% 
     group_by(target) %>%
     summarise(across(everything(), sum),
               .groups = 'drop') %>%
     as.data.frame() -> FPKtab_group
   
   
   smps <- c("L-CN.Prd", "L-RA.Prd", "M-CN.Prd", "M-RA.Prd",  "P-CN.Prd",
             "P-RA.Prd", "SA-RA.Prd", "SA-CN.Prd")  
   
   idx <- which(colnames(FPKtab_group) %in% smps)
   
   FPK_df <- FPKtab_group[,c(1,idx)]
   cnms <- gsub(".Prd", "", colnames(FPK_df))
   colnames(FPK_df) <- cnms
   rownames(FPK_df) <- FPK_df$target
   C_FPKM <-  FPK_df[,2:9]
   
   meta$Muestras <- colnames(C_FPKM)
   
   fit_data <- Maaslin2(
     C_FPKM, meta1, 'Out_Cazymes_maaslin2',
     #min_abundance = 0.005,
     fixed_effects = c('Treatment'), #P, Arena, Arcilla, Limo, CIC
     #random_effects = c('Sitio'), #'CO'
     plot_heatmap = TRUE,
     plot_scatter = TRUE,
     # standardize = FALSE,
     normalization = 'TSS'
   )
   
   
   ctb <- read_delim(file="Out_Cazymes_maaslin2/all_results.tsv")
   
   ctb %>%
     mutate(FC=exp(coef),
            log2FC=log2(FC),
            log10FC=log10(FC),
            log10qval=-log10(qval),
            log10pval=-log10(pval)) -> ctb1 
   
   ctb1$Land_use <- rep("NO", dim(ctb1)[1])
   ctb1$Land_use[which(ctb1$qval<0.25 & ctb1$coef<0)] <-"NG"
   ctb1$Land_use[which(ctb1$qval<0.25 & ctb1$coef>0)] <-"AR"
   
   ctb1$dalabel <- NA
   ctb1$dalabel[ctb1$Land_use != "NO"] <- ctb1$feature[ctb1$Land_use != "NO"]
   
   ctb1$dalabel <- gsub(".hmm", "", ctb1$dalabel)
   
   
   
   ggplot(data=ctb1, aes(x=log2FC, y=log10qval, label=dalabel, color=Land_use)) + 
     geom_point(size=2.5, alpha=.8) + 
     xlim(c(-10,10))+
     ylim(c(-0.5, 2.5))+
     theme_minimal()+
     theme(legend.position ='none')+
     geom_vline(xintercept=c(0), col="gray10", size=.4, alpha=.5) +
     geom_hline(yintercept=-log10(0.25), col="gray10", size=.4, alpha=.5)+
     geom_text_repel(size=2.8, max.overlaps = 10)+
     ggtitle('B')+
     ylab("log10(qval)")+
     xlab("log2(FoldChange)")+
     scale_color_manual(values=c('dodgerblue4','firebrick3' ,  "gray40")) -> Caz_vplot
   Caz_vplot
   

   ctb2 <- ctb1[which(!ctb1$Land_use=="NO"),]
   ctb2 <- ctb2[order(ctb2$Land_use),]
   sel <- c(1,4,5,8,9,10,11,13,15,16)
   ctb3 <- ctb2[,sel]
   
   ghx <- grep("GH",ctb3$feature)
   plx <- grep("PL",ctb3$feature)
   cex <- grep("CE",ctb3$feature)
   gtx <- grep("GT",ctb3$feature)
   cbx <- grep("CBM",ctb3$feature)
   
   class <- rep("Glycoside hydrolase", dim(ctb3)[1])
   class[plx] <- "Polysaccharide lyase"
   class[cex] <- "Carbohydrate esterase"
   class[gtx] <- "Glycosil transferase"
   class[cbx] <- "Carbohydrate binding module"
   
   ctb3$'Enzyme class' <- class
  
   #write_delim(ctb3, "CAZymes_supp_table.csv", delim=";")
  
   ###############################
   ### CAZymes substrate mapping
   ###############################
   
   tbs <- read_xls("Table S2.xls")
   colnames(tbs) <- tbs[1,]
   tbs <- tbs[2:744,]
   
   tfam <- tbs$Family
   tsub <- tbs$Substrate_high_level
   tec  <- tbs$EC_Number
   
   sbst <- character(0)
   
   for(j in 1:length(ctb3$feature)){
     
     caz <- gsub(".hmm", "", ctb3$feature[j])
     
     #remove subfamily
     czz <- strsplit(caz, "_")[[1]][1]
     
     #Compare with table
     idx <- which(tfam %in% caz)
     substrate <- unique(tsub[idx])
     
     #collapse substrates
     substrate2 <- paste(substrate, collapse="_")
     sbst <- c(sbst, substrate2)
     
   }
   
   
   ctb3$Substrate <- sbst
   
   coso <- separate_rows(ctb3,Substrate,sep="_")
   coso$Substrate <- gsub(",", "", coso$Substrate)
   
   sbst1 <- gsub("abscisic acid", "abscisic_acid", coso$Substrate)
   sbst2 <- gsub("arabinogalactan protein", "arabinogalactan_protein", sbst1)
   sbst3 <- gsub("host glycan", "host_glycan", sbst2)
   sbst4 <- gsub("sialic acid", "sialic_acid", sbst3)
   sbst5 <- gsub("cephalosporin C", "cephalosporin_C", sbst4)
   sbst6 <- gsub("human milk polysaccharide", " human_milk_polysaccharide", sbst5)
   sbst7 <- gsub("uric acid", "uric_acid", sbst6)
   
   coso$Substrate <- sbst7
   
   Table_final <- separate_rows(coso,Substrate,sep=" ")
   
   tbl_AR <- Table_final[which(Table_final$coef>0),]
   tbl_AR1 <- tbl_AR[which(!tbl_AR$Substrate==""),]
   tbl_NG <- Table_final[which(Table_final$coef<0),]
   tbl_NG1 <- tbl_NG[which(!tbl_NG$Substrate==""),]
   
  out<- c("ulvan", "porphyran", "palatinose", "glycosidases", "glucosylglycerol",
     "glucosylglycerate", "beta-glycan", "beta-glucuronan", "host_glycan", "", "agarose",
     "carrageenan")
   
   tbla <- Table_final[which(!Table_final$Substrate%in%out),]
   
   ggplot(data=tbla)+
     geom_bar(mapping=aes(x=forcats::fct_infreq(Substrate), fill=Land_use), stat="count", alpha=.85)+
       theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1, size=9),
             panel.background = element_rect(fill = "white", colour = "gray"),
             legend.position = "left")+
     #coord_flip()+
     xlab("Substrate")+
     ylab("Enzyme count")+
     ggtitle(label="C")+
     scale_fill_manual(values=c('dodgerblue4','firebrick3')) -> pcaz
    pcaz
    
    ######################################
    ## CAZymes Alpha diversity analysis ##
    ######################################
    
    phylocaz <- phyloseq(otu_table(C_FPKM, taxa_are_rows = TRUE), sample_data(meta1))
   
    df_rich <- as.data.frame(estimate_richness(phylocaz, split=TRUE, measures = "shannon"))
   
    trt  <- gsub("CN", "NG", meta1$Treatment)
    trt1 <- gsub("RA", "AR", trt)
    df_rich$Land_use <- as.factor(trt1)
    
    Sit <- gsub("Laurnaga", "L", meta1$Sitio)
    Siti <- gsub("Martirena", "M", Sit)
    Sitio <- gsub("Salto", "SA", Siti)
    Sitio1 <- gsub("Poni", "P", Sitio)
    
    meta1$Sitio <- Sitio1
    
    df_rich$Site <- meta1$Sitio
     
    alf_fun <- ggplot(data=df_rich, mapping = aes(x=as.factor(Land_use), y = Shannon, color = Land_use))+
      geom_jitter(mapping=aes(shape=Site), size=5, alpha=.6)+
      geom_boxplot(outlier.shape = NA, alpha = 0.1)+
      geom_signif(comparisons=list(c("NG", "AR")), map_signif_level = TRUE)+
      theme_classic()+
     # theme(legend.position = 'none')+
      xlab("Land use")+
      ylab("Shannon index")+
     # ylim(c(750, 1000))+
      ggtitle(label = "A")+
      scale_color_aaas()
    
    alf_fun

     ########################
    ###Microbial necromass recycling enzymes
   ########################
    
    chit <- Table_final[which(Table_final$Substrate=="chitin"),1]
    pept <- Table_final[which(Table_final$Substrate=="peptidoglycan"),1]
    
    chit_FPKM <- colSums(C_FPKM[which(rownames(C_FPKM)%in%chit$feature),])
    pept_FPKM <- colSums(C_FPKM[which(rownames(C_FPKM)%in%pept$feature),])
    
    Subs_FPKM <- as.data.frame(rbind(chit_FPKM, pept_FPKM))
    
    
    tab_FPKM <- pivot_longer(Subs_FPKM, cols = c("L-CN", "L-RA", "M-CN", "M-RA", "P-CN", "P-RA", "SA-CN", "SA-RA"))
    
    sit <- gsub("-CN", "", tab_FPKM$name )
    sit1 <- gsub("-RA", "", sit)
    
    cbind(tab_FPKM, 
          c(rep("Chitin", 8),rep("Petidoglycan",8)),
          rep(c("NG", "AR"), 8),
          sit1) -> ctab
      
   
    colnames(ctab) <- c("Samples", "FPKM", "Substrate", "Land_use", "Site")
   
    
    enz_fun <- ggplot(data=ctab, mapping = aes(x=Land_use, y = FPKM, color = Land_use))+
      geom_jitter(mapping=aes(shape=Site), size=4, alpha=.6)+
      geom_boxplot(outlier.shape = NA, alpha=.1)+
      geom_signif(comparisons=list(c("AR", "NG")), map_signif_level = TRUE, textsize=5)+
      facet_grid(~Substrate)+
      theme_classic()+
      xlab("Land use")+
      ylab("CAZYmes FPKM")+
      ggtitle(label = "C")+
      ylim(c(NA, 140))+
      scale_color_manual(values=c("#3A488AFF","#E5AD4FFF"))
    
    enz_fun

    
   ######################
   ###Final Figure layout 
   ######################
   
   gr1 <- plot_grid(ncyc_plot, pcyc_plot, ncol=2, rel_heights = c(1,1))
   
   gr2 <- plot_grid(alf_fun,Caz_vplot, pcaz,  ncol = 3, rel_widths = c(1,1,2), rel_heights = c(0.7,1,1)) 
   
   tiff("Fig_3.tiff", units="in", width=12, height=8, res=800)
  
   plot_grid(gr2, gr1, nrow = 2)  -> Fig3 
   Fig3  
   dev.off()
   
   


 
   
