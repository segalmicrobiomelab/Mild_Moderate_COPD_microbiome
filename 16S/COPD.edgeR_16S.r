#Load Packages
library(DESeq2)
library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(pathfindR)
library(scales)
library(data.table)
library(fBasics)
library(forcats)
library(omu)
library(maptools)
library(phyloseq)
library(vegan)
#library(ggpmisc)
library(dplyr)
library(tibble)
library(formattable)
library("htmltools")
library("webshot")    
library(splitstackshape)
library(decontam)
library(dplyr)
library(grid)
library(cowplot)
library(wesanderson)
library(colorspace)
library(factoextra)
library(Rtsne)
library(cluster)
library(umap)
library(magrittr)
library(ggpubr)
#library(egg)


#Set Theme for Figures
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")

#Choose Alpha/FDR
alpha = 0.2

#Create Function For Analysis
deseq <- function(metadata,counts,recounts,sampletype,comparison) {
        #Load MetaData
        #load(file="COPD.16S.RData")
        #Subset BAL BKG UA
        #BAL.BKG.UA.OTU.Table = subset_samples(OTU.Table, Sample_type_IS  %in% c("UA","BKG","BAL"))
        #BAL.BKG.UA.OTU.Rel.Table = subset_samples(OTU.Rel.Table, Sample_type_IS  %in% c("UA","BKG","BAL"))
        #Load MetaData
        coldata <- read.delim2(metadata, sep="\t")
        subjects <- read.delim2("NA_crew_forImran.txt", sep="\t")
        subjects2 <- data.frame(subjects$study_id,subjects$study_group)
        names(subjects2)<-c("Suject_ID2","Subject_Type_2")
        coldata <- merge(coldata,subjects2,by="Suject_ID2",all.x=TRUE)
        #coldata <- read.delim2("COPD_MAP_IS_2.txt", sep="\t")
        #Select only the samples with 16S data
        coldata <- coldata[coldata$sixteenS==1,]
        #replace all NAs
        #coldata$Subject_Type_2 <- ifelse(coldata$CAT_Score_High==1 & coldata$Subject_Type=="Smoker.Control","Symptomatic_SC",
        #                          ifelse(coldata$CAT_Score_High==0 & coldata$Subject_Type=="Smoker.Control","Asymptomatic_SC",
        #                          ifelse(coldata$CAT_Score_High=="n.a" & coldata$Subject_Type=="Smoker.Control","Asymptomatic_SC",
        #                          ifelse(coldata$Subject_Type=="COPD","COPD", NA))))
        coldata[ coldata == "n.a." ] <- NA
        coldata[ coldata == "n.a" ] <- NA
        #Create SampleID
        #coldata$SampleID <- coldata$SampleID_16S
        #coldata <- coldata %>% relocate(Description, .after = last_col())
        #coldata <- coldata%>%select(SampleID, everything())
        #coldata <- coldata[order(coldata$SampleID_16S),]       

        #Function to put NA in any empty cells
        empty_as_na <- function(x){
            if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
            ifelse(as.character(x)!="", x, NA)
            }
        ## transform all columns
        coldata <- coldata %>% mutate_each(funs(empty_as_na)) 
        coldata$Study_Linked_ID <- coldata$SampleID_16S
        coldata$Sample.Type <- coldata$Sample_type_IS
        coldata <- coldata[!duplicated(coldata$Study_Linked_ID), ]
        coldata <- coldata[coldata$Study_Linked_ID!="n.a", ]
        rownames(coldata) <- coldata$Study_Linked_ID
        #coldata$Subject_Type_2 <- ifelse(coldata$CAT_Score_High==1 & coldata$Subject_Type=="Smoker.Control","Symptomatic_SC",
        #                          ifelse(coldata$CAT_Score_High==0 & coldata$Subject_Type=="Smoker.Control","Asymptomatic_SC",
        #                          ifelse(coldata$Subject_Type=="COPD","COPD", NA)))

        #Load Count Data
        mycounts  <-read.delim2(paste0(counts,".txt"), sep="\t", row.names=1)
        #mycounts <-read.delim2("16S_COPD_OTU_Table_BAL_BKG_UA.txt", sep="\t", row.names=1)
        #Create Relative Abundance Table
        rel <-as.data.frame(mycounts) #make counts dataframe
        relcounts <-
            rel %>% 
                rownames_to_column('gs') %>%
                group_by(gs) %>% 
                summarise_all(funs(sum)) %>%
                mutate_if(is.numeric, funs(./sum(.))) %>%
                column_to_rownames('gs')
        #Find matching sample ID for both datasets
        needed<-which(rownames(coldata) %in% colnames(mycounts))    
        needed2<-which(rownames(coldata) %in% colnames(relcounts))    
        #keep only matching IDs from count data
        coldata2<-coldata[needed,]
        #Order Meta Data by SampleId
        coldata2 <- coldata2[order(coldata2$Study_Linked_ID),]
        #keep only matching IDs from count data
        wanted<-which(colnames(mycounts) %in% rownames(coldata))
        wanted2<-which(colnames(relcounts) %in% rownames(coldata))
        #keep only matching IDs from count data
        mycounts2<-mycounts[,wanted]
        relcounts2<-relcounts[,wanted2]
        #Order Count Data by SampleID
        mycounts2 <-mycounts2[, order(colnames(mycounts2))]
        relcounts2 <-relcounts2[, order(colnames(relcounts2))]
        #Convert any NAs to 0
        mycounts2[is.na(mycounts2)] <- 0
        relcounts2[is.na(relcounts2)] <- 0
        #Copy of Count Table
        mycounts3 <- mycounts2
        relcounts3 <- relcounts2
        #Convert Count Table into a Numeic Data Frame
        d1 = data.frame(lapply(mycounts3, function(x) as.numeric(as.character(x))),
                           check.names=F, row.names = rownames(mycounts3))
        d2 = data.frame(lapply(relcounts3, function(x) as.numeric(as.character(x))),
                           check.names=F, row.names = rownames(relcounts3))
        #Remove 0 counts in relative Table
        d2 <- d2[rowSums(d2[, -1] > 0) != 0, ]
        d2 <- d2 %>% select(which(!colSums(d2, na.rm=TRUE) %in% 0))
        #get the columns to match
        wanted<-which(colnames(d1) %in% colnames(d2))
        wanted2<-which(rownames(coldata2) %in% colnames(d2))
        d1<-d1[,wanted]
        coldata2<-coldata2[wanted2,]

        #Convert Data to Integers to Run DESEq
        d1[] <- lapply(d1, as.integer)
        count_sixteen <- d1
        meta_sixteen <- coldata2
        #Fix Rownames
        colnames(count_sixteen) <-gsub("\\.171","",colnames(count_sixteen))
        colnames(count_sixteen) <-gsub("\\.172","",colnames(count_sixteen))
        colnames(count_sixteen) <-gsub("\\.","_",colnames(count_sixteen))
        colnames(count_sixteen) <-gsub("BALF","BAL",colnames(count_sixteen))
        #Fix MetaData Colnames
        rownames(meta_sixteen) <-gsub("\\.171","",rownames(meta_sixteen))
        rownames(meta_sixteen) <-gsub("\\.172","",rownames(meta_sixteen))
        rownames(meta_sixteen) <-gsub("\\.","_",rownames(meta_sixteen))
        rownames(meta_sixteen) <-gsub("BALF","BAL",rownames(meta_sixteen))


        #Set Second DataSet Working Directory
        setwd("/gpfs/data/segallab/COPD/Metatranscriptome")
        coldata <- read.delim2("counts_taxa.metatranscriptome2.TSNE_Clustering.txt", sep="\t")
        #coldata <- read.delim2("KO_Pathway_broken_COPD_Metatranscriptome_sum_path03.TSNE_Clustering.txt",sep="\t")
        mycounts <-read.delim2("counts_taxa.metatranscriptome2.txt", sep="\t", row.names=1)
        #coldata <- read.delim2(metadata, sep="\t")
        coldata$Study_Linked_ID <- coldata$GTC_Updated_Name_MetatranscriptomeMH2021 #Metatranscriptome
        coldata$Study_Linked_ID <- gsub("-","_",coldata$Study_Linked_ID) #Metatranscriptome
        coldata$Study_Linked_ID <- gsub(".Sup","_Sup",coldata$Study_Linked_ID) #Metatranscriptome
        coldata <- coldata[!duplicated(coldata$Study_Linked_ID), ]
        coldata <- coldata[coldata$Study_Linked_ID!="n.a", ]
        rownames(coldata) <- coldata$Study_Linked_ID
        coldata$Sample.Type <- coldata$Sample_type_IS
        coldata <- coldata[!(coldata$Study_Linked_ID=="COPD_0002_BAL_L" | coldata$Study_Linked_ID=="COPD_0035_BAL_L" | coldata$Study_Linked_ID=="COPD_0041_BAL_L"),]
        #Load Count Data
        #mycounts  <-read.delim2(paste0(counts,".txt"), sep="\t", row.names=1)
        colnames(mycounts) <- gsub("\\.","_",colnames(mycounts))
        #Create Relative Abundance Table
        rel <-as.data.frame(mycounts) #make counts dataframe
        relcounts <-
            rel %>% 
                rownames_to_column('gs') %>%
                group_by(gs) %>% 
                summarise_all(funs(sum)) %>%
                mutate_if(is.numeric, funs(./sum(.))) %>%
                column_to_rownames('gs')
        #Find matching sample ID for both datasets
        needed<-which(rownames(coldata) %in% colnames(mycounts))    
        needed2<-which(rownames(coldata) %in% colnames(relcounts))    
        #keep only matching IDs from count data
        coldata2<-coldata[needed,]
        #Order Meta Data by SampleId
        coldata2 <- coldata2[order(coldata2$Study_Linked_ID),]
        #keep only matching IDs from count data
        wanted<-which(colnames(mycounts) %in% rownames(coldata))
        wanted2<-which(colnames(relcounts) %in% rownames(coldata))
        #keep only matching IDs from count data
        mycounts2<-mycounts[,wanted]
        relcounts2<-relcounts[,wanted2]
        #Order Count Data by SampleID
        mycounts2 <-mycounts2[, order(colnames(mycounts2))]
        relcounts2 <-relcounts2[, order(colnames(relcounts2))]
        #Convert any NAs to 0
        mycounts2[is.na(mycounts2)] <- 0
        relcounts2[is.na(relcounts2)] <- 0
        #Copy of Count Table
        mycounts3 <- mycounts2
        relcounts3 <- relcounts2
        #Convert Count Table into a Numeic Data Frame
        d1 = data.frame(lapply(mycounts3, function(x) as.numeric(as.character(x))),
                           check.names=F, row.names = rownames(mycounts3))
        d2 = data.frame(lapply(relcounts3, function(x) as.numeric(as.character(x))),
                           check.names=F, row.names = rownames(relcounts3))
        #Remove 0 counts in relative Table
        d2 <- d2[rowSums(d2[, -1] > 0) != 0, ]
        d2 <- d2 %>% select(which(!colSums(d2, na.rm=TRUE) %in% 0))
        #get the columns to match
        wanted<-which(colnames(d1) %in% colnames(d2))
        wanted2<-which(rownames(coldata2) %in% colnames(d2))
        d1<-d1[,wanted]
        coldata2<-coldata2[wanted2,]

        #Convert Data to Integers to Run DESEq
        d1[] <- lapply(d1, as.integer)
        count_metat <- d1
        meta_meta2 <- coldata2

        #keep only matching IDs from count data
        wanted<-which(colnames(count_sixteen) %in% colnames(count_metat))
        #keep only matching IDs from count data
        d1<-count_sixteen[,wanted]
        #Find matching sample ID for both datasets
        needed<-which(rownames(meta_sixteen) %in% colnames(d1))    
        #keep only matching IDs from count data
        coldata2<-meta_sixteen[needed,]
        rely <-
            d1 %>% 
                rownames_to_column('gs') %>%
                group_by(gs) %>% 
                summarise_all(funs(sum)) %>%
                mutate_if(is.numeric, funs(./sum(.))) %>%
                column_to_rownames('gs')
        idx <- rowSums(rely >= 0.00000000004 ) >= 2 #16S OTU
        rely_trim <- rely[idx,] #16S OTU
        #keep only matching IDs from count data
        wanted2<-which(rownames(d1) %in% rownames(rely_trim))
        #keep only matching IDs from count data
        #d1<-d1[wanted2,]
        ddsv <- DESeqDataSetFromMatrix(countData = d1,
                              colData = coldata2,
                              design= ~ Sample.Type)
        #Remove 0 counts from raw table                              
        #idx <- colSums( counts(ddsv)==0) 
        #ddsv <- ddsv[ , idx]
        #Remove non consented
        #ddsv <- ddsv[, ddsv$Paper_Final==1]                            
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsv), 1, gm_mean)
        #Estimate Factors of DESeq Object
        ddsv <- estimateSizeFactors(ddsv, geoMeans = geoMeans)
        ddsv1 <- ddsv[, ddsv$Sample.Type %in% c("BAL","UA")]
        ddsv2 <- ddsv[, ddsv$Sample.Type %in% c("BAL","BKG")]
        ddsv3 <- ddsv[, ddsv$Sample.Type %in% c("UA","BKG")]

        #Transforming data - Option #2 is variance stabilizing transformation
        vsdv <- varianceStabilizingTransformation(ddsv)
        vsdv1 <- vsdv[, vsdv$Sample.Type %in% c("BAL","UA")]
        vsdv2 <- vsdv[, vsdv$Sample.Type %in% c("BAL","BKG")]
        vsdv3 <- vsdv[, vsdv$Sample.Type %in% c("UA","BKG")]

        #vsdv <- rlog(ddsv, fitType="local",blind=FALSE)
        #----------------------
        ##PCOA PLOT
        #----------------------
        #Create Distance Matrix
        vsdv0 <- ifelse(assay(vsdv)<0,0,assay(vsdv))
        vsdv01 <- ifelse(assay(vsdv1)<0,0,assay(vsdv))
        vsdv02 <- ifelse(assay(vsdv2)<0,0,assay(vsdv))
        vsdv03 <- ifelse(assay(vsdv3)<0,0,assay(vsdv))

        vegdist   = vegdist(t(assay(ddsv)), method="bray")
        vegdist1   = vegdist(t(assay(ddsv1)), method="bray")
        vegdist2   = vegdist(t(assay(ddsv2)), method="bray")
        vegdist3   = vegdist(t(assay(ddsv3)), method="bray")

        #vsdv0 <- ifelse(d2<0,0,d2)
        #vegdist   = vegdist(t(d2), method="bray")
        #Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
        CmdScale <- cmdscale(vegdist, k =10)
        #calculated Sample variance for each PC
        vars <- apply(CmdScale, 2, var)
        #Create Variable with the Percent Variance
        percentVar <- round(100 * (vars/sum(vars)))
        #Merge PC Data with MetaData
        require(data.table)
        newResults <- merge(x = CmdScale, y = colData(vsdv), by = "row.names", all.x = TRUE)
        #Rename Variables for PC1 and PC2
        colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
        colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
        colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
        #Calculate the Centroid Value
        centroids <- aggregate(cbind(PC1,PC2)~ Sample.Type,data= newResults, mean)
        #Merge the Centroid Data into the PCOA Data
        newResults <- merge(newResults,centroids,by="Sample.Type",suffixes=c("",".centroid"))

        #Create Data For Statistics
        data.adonis <- data.frame(colData(vsdv))
        data.adonis1 <- data.adonis[data.adonis$Sample.Type %in% c("BAL","UA"),]
        data.adonis2 <- data.adonis[data.adonis$Sample.Type %in% c("BAL","BKG"),]
        data.adonis3 <- data.adonis[data.adonis$Sample.Type %in% c("UA","BKG"),]

        #Run the Statistics
        samplepermanova <- adonis(vegdist ~ Sample.Type, data.adonis)
        samplepermanova <- as.data.frame(samplepermanova$aov.tab)
        samplepermanova <- samplepermanova$'Pr(>F)'[1]

        samplepermanova1 <- adonis(vegdist1 ~ Sample.Type, data.adonis1)
        samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
        samplepermanova1 <- samplepermanova1$'Pr(>F)'[1] #0.001

        samplepermanova2 <- adonis(vegdist2 ~ Sample.Type, data.adonis2)
        samplepermanova2 <- as.data.frame(samplepermanova2$aov.tab)
        samplepermanova2 <- samplepermanova2$'Pr(>F)'[1] #0.001

        samplepermanova3 <- adonis(vegdist3 ~ Sample.Type, data.adonis3)
        samplepermanova3 <- as.data.frame(samplepermanova3$aov.tab)
        samplepermanova3 <- samplepermanova3$'Pr(>F)'[1] #0.001

        #PLOT IT
            ggsave(filename=paste0(counts,".BRAY_vsd_sampletype_PERMANOVA_",samplepermanova,".pdf"),
            ggplot(newResults, aes(PC1, PC2,color=Sample.Type)) + # Graph PC1 and PC2
            xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
            ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
            scale_color_manual(values=c("#1700F5", "#BEBEBE", "#F7A501")) + 
            geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Sample.Type))+ 
            geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=Sample.Type), size=10) + #BKG DFW Lung.Tissue.In Lung.Tissue.UnIn Mock
            geom_point(data=newResults,aes(color=Sample.Type),size=5,alpha=0.5) + # Set the size of the points
            theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
            panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
            panel.grid.minor = element_blank(),strip.background=element_blank(),
            axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
            axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),
            plot.margin=unit(c(1,1,1,1),"line"),legend.position="none"),
            height = 10, width = 10)
        #----------------------
        ##Cluster
        #----------------------
        #filter out genes where there are less than 3 samples with normalized counts greater than or equal to 100.
        #idx <- rowSums(relcounts3 >= 0.00004 ) >= 15 #Metatranscriptome All
        #df <- relcounts3[idx,] #Metatranscriptome All
        #idx <- rowSums(relcounts3 >= 0.0000004 ) >= 4 #16S OTU
        #df <- relcounts3[idx,] #16S OTU
        ##cluster Genuses(row)
        #umap_model_1 = umap(as.matrix(t(df)))
        #umapcord <- as.data.frame(umap_model_1$layout)
        #umapcord %<>% mutate(PC1 = umap_model_1$layout[, 1], PC2 = umap_model_1$layout[,2])
        #rownames(umapcord) <- rownames(coldata2)
        #tsne_model_1 = Rtsne(as.matrix(t(df)), check_duplicates=FALSE, pca=TRUE,theta=0.5,perplexity=5, dims=2)
        #tsnecord <- as.data.frame(tsne_model_1$Y)
        #tsnecord %<>% mutate(PC1 = tsne_model_1$Y[, 1], PC2 = tsne_model_1$Y[,2])
        #rownames(tsnecord) <- rownames(coldata2)
        ##Clustering
        #km.norm.tsne = kmeans(tsne_model_1$Y, 3, nstart = 25)
        #km.norm.umap = kmeans(umap_model_1$layout, 3, nstart = 25)
        #hc.norm.tsne = hclust(dist(tsne_model_1$Y)) #Complete
        #hc.norm.umap = hclust(dist(umap_model_1$layout)) #Complete
        ##Put the cluster classification in the PC table
        #tsnecord$kmeans = factor(km.norm.tsne$cluster)
        #umapcord$kmeans = factor(km.norm.umap$cluster)
        #tsnecord$hclust = factor(cutree(hc.norm.tsne, 3))
        #umapcord$hclust = factor(cutree(hc.norm.umap, 3))
        ##Merge TSNE with MetaData
        #    tsne <- merge(x = tsnecord, y = coldata2, by = "row.names", all.x = TRUE)
        #    tsne$col <- ifelse(tsne$kmeans==1,"Cluster_3", 
        #               ifelse(tsne$kmeans==2,"Cluster_2","Cluster_1"))
        #    tsne$col2 <- ifelse(tsne$hclust==1,"Cluster_3", 
        #               ifelse(tsne$hclust==2,"Cluster_2","Cluster_1"))            
        #    #Calculate the Centroid Value
        #    centroids <- aggregate(cbind(PC1,PC2)~ Sample.Type,data= tsne, mean)
        #    #Merge the Centroid Data into the PCOA Data
        #    tsne <- merge(tsne,centroids,by="Sample.Type",suffixes=c("",".centroid"))
        #    ggsave(filename=paste0(counts,".TSNE_Kmeans.PCOA.pdf"),
        #        ggplot(tsne, aes(PC1, PC2,color=Sample.Type)) + # Graph PC1 and PC2
        #        scale_color_manual(values=c("#1700F5", "#BEBEBE", "#F7A501")) + 
        #        geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Sample.Type))+ 
        #        geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=Sample.Type), size=10) + #BKG DFW Lung.Tissue.In Lung.Tissue.UnIn Mock
        #        geom_point(data=tsne,aes(fill=col),color="black",size=6, pch=21,alpha=0.5) + # Set the size of the points
        #        scale_fill_manual(values=c("#4DAC26","#FFD479","#D01C8B")) + 
        #        xlab("tsne 1")+
        #        ylab("tsne 2")+
        #        theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        #        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        #        panel.grid.minor = element_blank(),strip.background=element_blank(),
        #        axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        #        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),
        #        plot.margin=unit(c(1,1,1,1),"line"),legend.position="none"),
        #    height = 10, width = 10)
        #write.table(tsne,file=paste0(counts,".TSNE_Clustering.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE) 
        ##mapping.table=sample_data(read.table("MAP.txt", header=T, sep="\t", row.names=1))
        ##Create TSNE DataFrame
        #tsne <- as.data.frame(tsne)
        #rownames(tsne) <- tsne$Row.names
        ##Create Mapping Table
        #mapping.table=sample_data(tsne)
        ##Create OTUTable
        #OTUtable = otu_table(d1, taxa_are_rows=TRUE)
        ##Create Phyloseq Object
        #ps=phyloseq(otu_table(OTUtable), mapping.table)
        ###If you want to nomalize OTU table before
        ### To normalize data you need to set a function
        #normalizeSample = function(x) {
        #    x/sum(x)
        #        }
        #rel.ps = transformSampleCounts(ps, normalizeSample)
        ##Prune Data for BIPLOT
        #OTU.Rel.wh1 = genefilter_sample(rel.ps, filterfun_sample(function(x) x > 0.005), A = 0.005 * nsamples(rel.ps))
        #OTU.Rel.table1B = prune_taxa(OTU.Rel.wh1, rel.ps)
        ##Ordinate NMDS and BRAY
        #OTU.ords <- ordinate(OTU.rel.table1B , "NMDS", "bray")
        ##Gather Data For Biplot from Ordination
        #pord <- plot_ordination(OTU.rel.table1B, OTU.ords, type = "biplot", color = "Genus", title = "Lower Airway Biplot", shape="Sample.Type")
        ##pord <- plot_ordination(OTU.rel.table1B, OTU.ords, type = "biplot", color = "Genus", title = "Lower Airway Biplot", shape = "Neutrophil_Elastase_Cat")
        ##Put Data into a dataframe
        #ord <- as.data.frame(pord$data)
        ##Select only the axes data
        #df <- subset(ord,select=c("NMDS1","NMDS2"))
        ##calculated Sample variance for each PC
        #vars <- apply(df, 2, var)
        ##Create Variable with the Percent Variance
        #percentVar <- round(100 * (vars/sum(vars)))
#
        ##Rename Variables for PC1 and PC2
        #colnames(ord)[colnames(ord)=="NMDS1"] <- "PC1"
        #colnames(ord)[colnames(ord)=="NMDS2"] <- "PC2"
#
        ##decide what otu to save 
        #otu.to.save <-as.character(ord$Strain)
        ##from relative table we should get the mean across the row of the otu table
        #OTU.rel.table.df <- data.frame(otu_table(x10))
        #OTU.rel.table.df.meanRA <- rowMeans(OTU.rel.table.df)
        ##need to subset AND reorder just the otus that we have 
        #OTU.rel.table.df.meanRA.save <- OTU.rel.table.df.meanRA[otu.to.save]
        ##add the abundnace data for the res dataframe
        #ord$abundance <- OTU.rel.table.df.meanRA.save
        ##Calculate the Centroid Value
        #centroids <- aggregate(cbind(PC1,PC2)~Sample.Type,data= ord, mean)
        ##Merge the Centroid Data into the PCOA Data
        #ord <- merge(ord,centroids,by="Sample.Type",suffixes=c("",".centroid"))
        ##Create Taxa name with Genus and OTU
        #ord$taxa <- paste(ord$Genus,ord$Strain,sep="_")
        #ord$taxa2 <- ifelse(ord$taxa=="Samples_NA",ord$taxa,as.character(ord$Genus))
        #pdf("16S_SampleType_Cluster_BiPlot.pdf", height = 15, width = 20)
        #    ggplot(ord, aes(PC1, PC2, color=Sample.Type)) +
        #    geom_point(size= ifelse(ord$Sample.Type=="Taxa" & ord$abundance>0.00037, 200 * ord$abundance, 0),alpha=0.7) +    
        #    geom_point(data=subset(ord,Sample.Type!="Taxa"),size=5,alpha=0.7) +
        #    #geom_point(data=subset(ord,FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa"),size= ifelse(FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa", 200 * ord$abundance, 2)) +
        #    xlab(paste0("NMDS1: ",percentVar[1],"% variance")) +
        #    ylab(paste0("NMDS2: ",percentVar[2],"% variance")) + 
        #    #coord_fixed() +
        #    #scale_color_manual(values=c("#BEBEBE","#FF5DE7", "#EA3323","#00CED1","#296218")) + 
        #    #plot ellipse
        #    #stat_ellipse(type = "t") + 
        #    #plot point and lines from centroid
        #    geom_point(data=subset(centroids,Sample.Type!="Taxa"), aes(x=PC1, y=PC2, color=Sample.Type), size=0) +
        #    geom_segment(data=subset(ord,Sample.Type!="Taxa"), aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Sample.Type))+ 
        #    #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
        #    #labels centroids 
        #    #geom_text_repel(aes(label=ifelse(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB=="Pseudomonas" , as.character(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB),'')),size=3,force=25) +
        #    geom_label_repel(data=subset(centroids,Sample.Type!="Taxa"), aes(x=PC1, y=PC2, label=c("BAL", "BKG", "UA")), size=10) +
        #    geom_text_repel(data=subset(ord,Sample.Type=="Taxa" & ord$abundance>0.00037),aes(label=taxa), size=4)+
        #    #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=lab), parse=TRUE,size=10) +
        #    #scale_x_reverse() +
        #    theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
        #dev.off()
        ##Subset Rel Table for BAL only
        #wanted<-which(rownames(d1) %in% rownames(df))
        ##keep only matching IDs from count data
        #d3<-d1[wanted,]
        #Create Deseq object for BAL analysis
        ddsvbal <- DESeqDataSetFromMatrix(countData = d1,
                              colData = coldata2,
                              design= ~ Sample.Type)
        #Subset just the BAL
        ddsvbal <- ddsvbal[, ddsvbal$Sample.Type %in% "BAL"]
        #ddsvbal <- ddsvbal[, ddsvbal$Subject_Type_2 %in% c("COPD","Asymptomatic_SC","Symptomatic_SC")]
        #Subset Rel Table for BAL only
        #wanted<-which(colnames(d2) %in% colnames(ddsvbal))
        #keep only matching IDs from count data
        #d3<-d2[,wanted]
        #Covert Variable to Factor
        ddsvbal$three_groups <- as.factor(ddsvbal$Subject_Type_2)
        #Create Seperate tables with each of the different compairsons
        #ddsvbal1 <- ddsvbal[, ddsvbal$three_groups %in% c("COPD","Asymptomatic_SC")]
        #ddsvbal2 <- ddsvbal[, ddsvbal$three_groups %in% c("Asymptomatic_SC","Symptomatic_SC")]
        #ddsvbal3 <- ddsvbal[, ddsvbal$three_groups %in% c("Symptomatic_SC","COPD")]
        #Create Tables of First Comparison
        balcounts <- assay(ddsvbal)
        balmeta   <- colData(ddsvbal)
        #Write Taxa Data
        #Get Phyloseq Data
        setwd("/gpfs/data/segallab/COPD/16S")
        load(file="COPD.16S.KMEANS.RData")
        #Subset BAL Table
        BAL.OTU.Table = subset_samples(OTU.Table, Sample_type_IS  %in% c("BAL"))
        #OTU.Rel.Table = transformSampleCounts(OTU.Table, normalizeSample)
        #BAL.OTU.Rel.Table = subset_samples(OTU.Rel.Table, Sample_type_IS  %in% c("BAL"))
        #Get Taxa Names from Phyloseq Object
        merged = cbind(as(balcounts, "data.frame"), as(tax_table(BAL.OTU.Table)[rownames(balcounts), ], "matrix"))
        #Convert Resuts table into a data.frame
        merged <- as.data.frame(merged)
        #convert to character
        merged$ASV <- rownames(merged)
        #Replace any no genus annotation as NA
        merged[merged==" s__"]<-NA
        merged[merged==" g__"]<-NA
        merged[merged==" f__"]<-NA
        merged[merged==" o__"]<-NA
        merged[merged==" c__"]<-NA
        merged[merged=="s__"]<-NA
        merged[merged=="g__"]<-NA
        merged[merged=="f__"]<-NA
        merged[merged=="o__"]<-NA
        merged[merged=="c__"]<-NA
        #Create name with family and (u.g)
        merged$gs <- ifelse(is.na(merged$Species),paste0(merged$Genus,"_",merged$ASV), paste0(merged$Genus,"_",merged$Species,"_",merged$ASV))
        merged$gs <- ifelse(is.na(merged$Genus),paste0(merged$Family, "_",merged$ASV), merged$gs)
        merged$gs <- ifelse(is.na(merged$Family), paste0(merged$Order,"_",merged$ASV),merged$gs)
        merged$gs <- ifelse(is.na(merged$Order), paste0(merged$Class, "_",merged$ASV),merged$gs)
        merged$gs <- ifelse(is.na(merged$Class), paste0(merged$Phylum,"_",merged$ASV),merged$gs)
        #Make the full trail the First Column
        merged$names <- merged$ASV
        merged$Gene.symbol <- merged$gs 
        merged$NAME <- paste0(merged$Genus,"_",merged$Species)    
        
    
        #Drop unnecessary Variables
        drops <- c("Domain","Phylum","Class","Order","Family","Genus","OTU","gs","Species","ASV","Kingdom","names","row2")
        res <- res[ , !(names(res) %in% drops)]
        #write.table(balcounts,file="metagenome_bacteria_BAL.txt", sep="\t", col.names = NA, row.names = TRUE)
        ##Create Tables of First Comparison
        #balcounts1 <- assay(ddsvbal1)
        #balmeta1   <- colData(ddsvbal1)
        ###Create Tables of Second Comparison
        #balcounts2 <- assay(ddsvbal2)
        #balmeta2   <- colData(ddsvbal2)
        ###Create Tables of Third Comparison
        #balcounts3 <- assay(ddsvbal3)
        #balmeta3   <- colData(ddsvbal3)
        #Create new DESEQ2 Objects
        ddsvbal <- DESeqDataSetFromMatrix(countData = balcounts,
                                      colData = balmeta,
                                      design= ~ Subject_Type)
        #ddsvbal1 <- DESeqDataSetFromMatrix(countData = balcounts1,
        #                              colData = balmeta1,
        #                              design= ~ three_groups)
        #ddsvbal2 <- DESeqDataSetFromMatrix(countData = balcounts2,
        #                              colData = balmeta2,
        #                              design= ~ three_groups)
        #ddsvbal3 <- DESeqDataSetFromMatrix(countData = balcounts3,
        #                            colData = balmeta3,
        #                            design= ~ three_groups)
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsvbal), 1, gm_mean)        
        #Estimate Factors of DESeq Object
        ddsvbal <- estimateSizeFactors(ddsvbal, geoMeans = geoMeans)
        #Calculate geometric means prior to estimate size factor
        #gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        #geoMeans = apply(counts(ddsvbal1), 1, gm_mean)
        ###Estimate Factors of DESeq Object
        #ddsvbal1 <- estimateSizeFactors(ddsvbal1, geoMeans = geoMeans)
        ###Calculate geometric means prior to estimate size factor
        #gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        #geoMeans = apply(counts(ddsvbal2), 1, gm_mean)
        ###Estimate Factors of DESeq Object
        #ddsvbal2 <- estimateSizeFactors(ddsvbal2, geoMeans = geoMeans)
        ###Calculate geometric means prior to estimate size factor
        #gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        #geoMeans = apply(counts(ddsvbal3), 1, gm_mean)
        ###Estimate Factors of DESeq Object
        #ddsvbal3 <- estimateSizeFactors(ddsvbal3, geoMeans = geoMeans)
        #Variance Stablize the data
        vsdvbal <- varianceStabilizingTransformation(ddsvbal)
        #vsdvbal1 <- varianceStabilizingTransformation(ddsvbal1)
        #vsdvbal2 <- varianceStabilizingTransformation(ddsvbal2)
        #vsdvbal3 <- varianceStabilizingTransformation(ddsvbal3)
        #vsdvbal <- rlog(ddsvbal, fitType="local",blind=FALSE)
        #----------------------
        ##PCOA
        #----------------------
        #Create Distance Matrix
        vsdvbal0 <- ifelse(assay(vsdvbal)<0,0,assay(vsdvbal))
        vegdist   = vegdist(t(vsdvbal0), method="bray")
        #vegdist   = vegdist(t(assay(ddsvbal3)), method="bray")
        #vegdist   = vegdist(t(d3), method="bray")
        #Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
        CmdScale <- cmdscale(vegdist, k =10)
        #calculated Sample variance for each PC
        vars <- apply(CmdScale, 2, var)
        #Create Variable with the Percent Variance
        percentVar <- round(100 * (vars/sum(vars)))
        #Merge PC Data with MetaData
        require(data.table)
        newResults <- merge(x = CmdScale, y = colData(vsdvbal), by = "row.names", all.x = TRUE)
        #Rename Variables for PC1 and PC2
        colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
        colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
        colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
        #Calculate the Centroid Value
        centroids <- aggregate(cbind(PC1,PC2)~ three_groups,data= newResults, mean)
        #Merge the Centroid Data into the PCOA Data
        newResults <- merge(newResults,centroids,by="three_groups",suffixes=c("",".centroid"))
        #Create Table for Statistics    
        data.adonis <- data.frame(colData(vsdvbal))
        #data.adonis <- data.frame(colData(vsdvbal1))
        #data.adonis <- data.frame(colData(vsdvbal2))
        #data.adonis <- data.frame(colData(vsdvbal3))
        #Run the Statistics
        samplepermanova <- adonis(vegdist ~ three_groups, data.adonis)
        samplepermanova <- as.data.frame(samplepermanova$aov.tab)
        samplepermanova <- samplepermanova$'Pr(>F)'[1]
        #PLOT IT
            ggsave(filename=paste0(counts,".BRAY_vsd_Subject_Type_2_PERMANOVA_",samplepermanova,".pdf"),
            ggplot(newResults, aes(PC1, PC2,color=three_groups)) + # Graph PC1 and PC2
            xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
            ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
            scale_color_manual(values=c("Asymptomatic_SC"="#4DAC26","Symptomatic_SC"="#FFD479","COPD"="#D01C8B")) + 
            geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= three_groups))+ 
            geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=three_groups), size=10) + #BKG DFW Lung.Tissue.In Lung.Tissue.UnIn Mock
            geom_point(data=newResults,aes(color=three_groups),size=5,alpha=0.5) + # Set the size of the points
            theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
            panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
            panel.grid.minor = element_blank(),strip.background=element_blank(),
            axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
            axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),
            plot.margin=unit(c(1,1,1,1),"line"), legend.position="none"),
            height = 10, width = 10)

        #Loadings
        #Set Variables
        counts <- balcounts
        groups <- as.data.frame(balmeta)
        vsdvbal00 <- as.data.frame(vsdvbal0)
        load(file="COPD.16S.KMEANS.RData")
        #Subset BAL Table
        BAL.OTU.Table = subset_samples(OTU.Table, Sample_type_IS  %in% c("BAL"))
        #OTU.Rel.Table = transformSampleCounts(OTU.Table, normalizeSample)
        #BAL.OTU.Rel.Table = subset_samples(OTU.Rel.Table, Sample_type_IS  %in% c("BAL"))
        #Get Taxa Names from Phyloseq Object
        merged = cbind(as(vsdvbal00, "data.frame"), as(tax_table(BAL.OTU.Table)[rownames(vsdvbal00), ], "matrix"))
        #Convert Resuts table into a data.frame
        merged <- as.data.frame(merged)
        #convert to character
        merged$ASV <- rownames(merged)
        #Replace any no genus annotation as NA
        merged[merged==" s__"]<-NA
        merged[merged==" g__"]<-NA
        merged[merged==" f__"]<-NA
        merged[merged==" o__"]<-NA
        merged[merged==" c__"]<-NA
        merged[merged=="s__"]<-NA
        merged[merged=="g__"]<-NA
        merged[merged=="f__"]<-NA
        merged[merged=="o__"]<-NA
        merged[merged=="c__"]<-NA
        #Create name with family and (u.g)
        merged$gs <- ifelse(is.na(merged$Species),paste0(merged$Genus,"_",merged$ASV), paste0(merged$Genus,"_",merged$Species,"_",merged$ASV))
        merged$gs <- ifelse(is.na(merged$Genus),paste0(merged$Family, "_",merged$ASV), merged$gs)
        merged$gs <- ifelse(is.na(merged$Family), paste0(merged$Order,"_",merged$ASV),merged$gs)
        merged$gs <- ifelse(is.na(merged$Order), paste0(merged$Class, "_",merged$ASV),merged$gs)
        merged$gs <- ifelse(is.na(merged$Class), paste0(merged$Phylum,"_",merged$ASV),merged$gs)
        #Make the full trail the First Column
        merged$names <- merged$ASV
        rownames(merged) <- merged$gs     
        #Drop unnecessary Variables
        drops <- c("Domain","Phylum","Class","Order","Family","Genus","OTU","gs","Species","ASV","Kingdom","names","row2")
        merged <- merged[ , !(names(merged) %in% drops)]
        ###
        wascores <- wascores(CmdScale, as.data.frame(t(merged)))
        wascores.comp1 <- wascores[order(abs(wascores[,1]),decreasing=T),1][1:20]
        wascores.comp2 <- wascores[order(abs(wascores[,2]),decreasing=T),2][1:20]
           
        plot.wascores.comp1 <- ggplot()+
          geom_bar(aes(x=!!wascores.comp1,
                       y=!!factor(names(wascores.comp1),levels=rev(names(wascores.comp1)))),
                   stat="identity")+    
          xlab("weighted average scores")+
          ylab("Taxa")+
          theme

        plot.wascores.comp2 <- ggplot()+
          geom_bar(aes(x=!!wascores.comp2,
                       y=!!factor(names(wascores.comp2),levels=rev(names(wascores.comp2)))),
                   stat="identity")+
          xlab("weighted average scores")+
          ylab("Taxa")+
          theme

        library(cowplot)
        pdf("16S.BAL_PCOA_Loadings.pdf", height = 5, width = 25)
            plot_grid(plot.wascores.comp1, plot.wascores.comp2, labels = c('A', 'B'))    
        dev.off()
        ##----------------------
        ###Alpha Diversity
        ##----------------------
        ##Calcultes Shannon Diversity
        #ddsvbal$Shannon = diversity(vsdvbal0, index = "shannon", MARGIN = 2, base = exp(1))
        ##ddsvbal$Shannon = diversity(d3, index = "shannon", MARGIN = 2, base = exp(1))
        ##Convert to data frame for ggplot
        #shannon = as.data.frame(colData(ddsvbal))
        ##Set Order Of Figure
        #shannon$or <-ifelse(shannon$three_groups=="Asymptomatic_SC", 1,NA)
        #shannon$or <-ifelse(shannon$three_groups=="Symptomatic_SC", 2,shannon$or)
        #shannon$or <-ifelse(shannon$three_groups=="COPD",3 ,shannon$or)
        ##Make Sure Shannon is Numeric
        #shannon$Shannon <- as.numeric(as.character(shannon$Shannon))
        ##Make sure Group is factor
        #shannon$Subject_Type <- as.factor(shannon$Subject_Type_2)
        #shannon1 <- shannon[shannon$three_groups %in% c("Asymptomatic_SC","Symptomatic_SC"),]
        #kruskal.test(Shannon ~ three_groups, data = shannon1)
        ##Kruskal-Wallis chi-squared = 0.8127, df = 1, p-value = 0.5516
        #shannon2 <- shannon[shannon$three_groups %in% c("COPD","Symptomatic_SC"),]
        #kruskal.test(Shannon ~ three_groups, data = shannon2)
        ##Kruskal-Wallis chi-squared = 0.87491, df = 1, p-value = 0.879
        #shannon3 <- shannon[shannon$three_groups %in% c("Asymptomatic_SC","COPD"),]
        #kruskal.test(Shannon ~ three_groups, data = shannon3)
        ##Kruskal-Wallis chi-squared = 0.031559, df = 1, p-value = 0.5126
#
        ##Check Statistics
        #sampleshannon <- kruskal.test(Shannon ~ three_groups, data = shannon)
        #sampleshannon <- sampleshannon$p.value
        ##PLOT IT
        #    ggsave(filename=paste0(counts,".SHANNON_vsd_Subject_Type2_Pvalue_",sampleshannon,".pdf"),
        #    ggplot(shannon, aes(x= reorder(three_groups, +or), y=Shannon, fill=three_groups)) + 
        #    stat_boxplot(geom ='errorbar', width=0.1)+
        #    geom_boxplot(outlier.shape = NA, width=0.5)+
        #    geom_jitter(shape=1, position=position_jitter(0.2))+
        #    scale_fill_manual(values=c("Asymptomatic_SC"="#4DAC26","Symptomatic_SC"="#FFD479","COPD"="#D01C8B")) + 
        #    #scale_x_discrete(labels = c('<28 Days','>28 Days','Deceased'))+ 
        #    ylab("Shannon Diversity") + 
        #    xlab("Subject Type")+
        #    #scale_y_continuous(breaks = seq(2, 5, .5), limits = c(2, 5)) + #For Fungi
        #    #scale_y_continuous(breaks = seq(3, 5, .5), limits = c(3, 5)) + #For Fungi
        #    theme,
        #    useDingbats=FALSE,
        #    height = 7, width = 5)
        #----------------------
        ##Differential Analysis
        #----------------------
#DropLevels
        ddsvbal$Subject_Type <- droplevels(ddsvbal$Subject_Type)
        #ddsvbal$three_groups <- droplevels(ddsvbal$Subject_Type_2)
        #ddsvbal1$three_groups <- droplevels(ddsvbal1$three_groups)
        #ddsvbal2$three_groups <- droplevels(ddsvbal2$three_groups)
        #ddsvbal3$three_groups <- droplevels(ddsvbal3$three_groups)
        #Set Reference
        ddsvbal$Subject_Type <- relevel(ddsvbal$Subject_Type, ref ="Smoker.Control")
        #ddsvbal1$three_groups <- relevel(ddsvbal1$three_groups, ref ="Asymptomatic_SC")
        #ddsvbal2$three_groups <- relevel(ddsvbal2$three_groups, ref ="Asymptomatic_SC")
        #ddsvbal3$three_groups <- relevel(ddsvbal3$three_groups, ref ="Symptomatic_SC")
        #Run the differential Analysis: Lung Cancer Vs Wild Type --> positive is upregulated in Lung Cancer; Negative is down regulated
        ddsvbal  <- DESeq(ddsvbal)
        #ddsvbal1  <- DESeq(ddsvbal1)
        #ddsvbal2  <- DESeq(ddsvbal2)
        #ddsvbal3  <- DESeq(ddsvbal3)
        #Output Result Table
        res     <- results(ddsvbal, cooksCutoff=FALSE)
        #res1     <- results(ddsvbal1, cooksCutoff=FALSE)
        #res2     <- results(ddsvbal2, cooksCutoff=FALSE)
        #res3     <- results(ddsvbal3, cooksCutoff=FALSE)
        #----------------------
        ##TAXA DATA
        #----------------------
        coldata3 <- meta_sixteen
        #Get Phyloseq Data
        setwd("/gpfs/data/segallab/COPD/16S")
        load(file="COPD.16S.KMEANS.RData")
        #Subset BAL Table
        BAL.OTU.Table = subset_samples(OTU.Table, Sample_type_IS  %in% c("BAL"))
        #OTU.Rel.Table = transformSampleCounts(OTU.Table, normalizeSample)
        #BAL.OTU.Rel.Table = subset_samples(OTU.Rel.Table, Sample_type_IS  %in% c("BAL"))
        #Get Taxa Names from Phyloseq Object
        res = cbind(as(res, "data.frame"), as(tax_table(BAL.OTU.Table)[rownames(res), ], "matrix"))
        #res1 = cbind(as(res1, "data.frame"), as(tax_table(BAL.OTU.Table)[rownames(res1), ], "matrix"))
        #res2 = cbind(as(res2, "data.frame"), as(tax_table(BAL.OTU.Table)[rownames(res2), ], "matrix"))
        #res3 = cbind(as(res3, "data.frame"), as(tax_table(BAL.OTU.Table)[rownames(res3), ], "matrix"))
        ##Create contam table
        #res4 <- tax_table(BAL.OTU.Table)
        #res4 <- as.data.frame(res4)
        ##res4 <- otu_table(BAL.OTU.Rel.Table)
        ##res4 <- as.data.frame(res4)
        ##res4 = cbind(as(res4, "data.frame"), as(tax_table(BAL.OTU.Rel.Table)[rownames(res4), ], "matrix"))
        ##convert to character
        #res4$ASV <- rownames(res4)
        ##Replace any no genus annotation as NA
        #res4[res4==" s__"]<-NA
        #res4[res4==" g__"]<-NA
        #res4[res4==" f__"]<-NA
        #res4[res4==" o__"]<-NA
        #res4[res4==" c__"]<-NA
        #res4[res4=="s__"]<-NA
        #res4[res4=="g__"]<-NA
        #res4[res4=="f__"]<-NA
        #res4[res4=="o__"]<-NA
        #res4[res4=="c__"]<-NA
        ##Remove the g_ and s_
        ##res4$Genus <- gsub("g__","",res4$Genus,fixed = TRUE)   
        ##res4$Species <- gsub("s__","",res4$Species,fixed = TRUE)   
        ##res4$Species <- gsub("\\s","",res4$Species,fixed = TRUE)   
        ##Create name with family and (u.g)
        #res4$gs <- ifelse(is.na(res4$Species),paste0(res4$Genus,"_",res4$ASV), paste0(res4$Genus,"_",res4$Species,"_",res4$ASV))
        #res4$gs <- ifelse(is.na(res4$Genus),paste0(res4$Family, "_",res4$ASV), res4$gs)
        #res4$gs <- ifelse(is.na(res4$Family), paste0(res4$Order,"_",res4$ASV),res4$gs)
        #res4$gs <- ifelse(is.na(res4$Order), paste0(res4$Class, "_",res4$ASV),res4$gs)
        #res4$gs <- ifelse(is.na(res4$Class), paste0(res4$Phylum,"_",res4$ASV),res4$gs)
        ##Make the full trail the First Column
        #res4$taxa <- res4$ASV
        #res4$Gene.symbol <- res4$gs
        #res4$taxa <- paste0(res4$Genus,"_",res4$Species)
        #res4$taxa <- gsub("g__","",res4$taxa,fixed = TRUE)   
        #res4$taxa <- gsub("s__","",res4$taxa,fixed = TRUE)   
        #res4$taxa <- gsub("NA","",res4$taxa,fixed = TRUE)   
        #res4$taxa <- gsub(" ","",res4$taxa,fixed = TRUE)   
        #res4$taxa <- gsub("_"," ",res4$taxa,fixed = TRUE)
        #res4$gss <- paste0(res4$Genus,"_",res4$Species)
        #res5 <- res4 %>% select(taxa,Gene.symbol)
        #res5 <- res5 %>% mutate(across(c("taxa"), ~ifelse(.=="", NA, as.character(.))))        
        #Vector of Colors
        #Get contam Data
        #contam <- read.delim2("contaminant16S_ASV.txt", sep="\t")
        #contam2 <- data.frame(contam$taxa,contam$contaminant)
        #names(contam2) <- c("taxa","contaminant")
        ##contam2 <-  separate(contam2, taxa, into = c("Genus", "Species"), sep = " (?=[^ ]+$)")
        ##contam2$taxa <- paste0("g__",contam2$Genus,"_","s__",contam2$Species)
        ##Merge with DESEQ data
        #contam2 <- merge(res4,contam2,by="taxa",all.x=TRUE)
        #contam2 <- contam2 %>% select(Gene.symbol,contaminant)
        #Create Color vector for tick labels
        #colvec <- ifelse(resy$probable_metatranscriptome_contaminant=="TRUE", "red", "black")
        #Fix the Taxa Names
        #Convert Resuts table into a data.frame
        res <- as.data.frame(res)
        #convert to character
        res$ASV <- rownames(res)
        #Replace any no genus annotation as NA
        res[res==" s__"]<-NA
        res[res==" g__"]<-NA
        res[res==" f__"]<-NA
        res[res==" o__"]<-NA
        res[res==" c__"]<-NA
        res[res=="s__"]<-NA
        res[res=="g__"]<-NA
        res[res=="f__"]<-NA
        res[res=="o__"]<-NA
        res[res=="c__"]<-NA
        #Create name with family and (u.g)
        res$gs <- ifelse(is.na(res$Species),paste0(res$Genus,"_",res$ASV), paste0(res$Genus,"_",res$Species,"_",res$ASV))
        res$gs <- ifelse(is.na(res$Genus),paste0(res$Family, "_",res$ASV), res$gs)
        res$gs <- ifelse(is.na(res$Family), paste0(res$Order,"_",res$ASV),res$gs)
        res$gs <- ifelse(is.na(res$Order), paste0(res$Class, "_",res$ASV),res$gs)
        res$gs <- ifelse(is.na(res$Class), paste0(res$Phylum,"_",res$ASV),res$gs)
        #Make the full trail the First Column
        res$names <- res$ASV
        res$Gene.symbol <- res$gs     
        
        #Convert Resuts table into a data.frame
        #res1 <- as.data.frame(res1)
        ##convert to character
        #res1$ASV <- rownames(res1)
        ##Replace any no genus annotation as NA
        #res1[res1==" s__"]<-NA
        #res1[res1==" g__"]<-NA
        #res1[res1==" f__"]<-NA
        #res1[res1==" o__"]<-NA
        #res1[res1==" c__"]<-NA
        #res1[res1=="s__"]<-NA
        #res1[res1=="g__"]<-NA
        #res1[res1=="f__"]<-NA
        #res1[res1=="o__"]<-NA
        #res1[res1=="c__"]<-NA
        ##Create name with family and (u.g)
        #res1$gs <- ifelse(is.na(res1$Species),paste0(res1$Genus,"_",res1$ASV), paste0(res1$Genus,"_",res1$Species,"_",res1$ASV))
        #res1$gs <- ifelse(is.na(res1$Genus),paste0(res1$Family, "_",res1$ASV), res1$gs)
        #res1$gs <- ifelse(is.na(res1$Family), paste0(res1$Order,"_",res1$ASV),res1$gs)
        #res1$gs <- ifelse(is.na(res1$Order), paste0(res1$Class, "_",res1$ASV),res1$gs)
        #res1$gs <- ifelse(is.na(res1$Class), paste0(res1$Phylum,"_",res1$ASV),res1$gs)
        ##Make the full trail the First Column
        #res1$names <- res1$ASV
        #res1$Gene.symbol <- res1$gs        
        #
        ##Convert Resuts table into a data.frame
        #res2 <- as.data.frame(res2)
        ##convert to character
        #res2$ASV <- rownames(res2)
        ##Replace any no genus annotation as NA
        #res2[res2==" s__"]<-NA
        #res2[res2==" g__"]<-NA
        #res2[res2==" f__"]<-NA
        #res2[res2==" o__"]<-NA
        #res2[res2==" c__"]<-NA
        #res2[res2=="s__"]<-NA
        #res2[res2=="g__"]<-NA
        #res2[res2=="f__"]<-NA
        #res2[res2=="o__"]<-NA
        #res2[res2=="c__"]<-NA
        ##Create name with family and (u.g)
        #res2$gs <- ifelse(is.na(res2$Species),paste0(res2$Genus,"_",res2$ASV), paste0(res2$Genus,"_",res2$Species,"_",res2$ASV))
        #res2$gs <- ifelse(is.na(res2$Genus),paste0(res2$Family, "_",res2$ASV), res2$gs)
        #res2$gs <- ifelse(is.na(res2$Family), paste0(res2$Order,"_",res2$ASV),res2$gs)
        #res2$gs <- ifelse(is.na(res2$Order), paste0(res2$Class, "_",res2$ASV),res2$gs)
        #res2$gs <- ifelse(is.na(res2$Class), paste0(res2$Phylum,"_",res2$ASV),res2$gs)
        ##Make the full trail the First Column
        #res2$names <- res2$ASV
        #res2$Gene.symbol <- res2$gs  
        #
        ##Convert Resuts table into a data.frame
        #res3 <- as.data.frame(res3)
        ##convert to character
        #res3$ASV <- rownames(res3)
        ##Replace any no genus annotation as NA
        #res3[res3==" s__"]<-NA
        #res3[res3==" g__"]<-NA
        #res3[res3==" f__"]<-NA
        #res3[res3==" o__"]<-NA
        #res3[res3==" c__"]<-NA
        #res3[res3=="s__"]<-NA
        #res3[res3=="g__"]<-NA
        #res3[res3=="f__"]<-NA
        #res3[res3=="o__"]<-NA
        #res3[res3=="c__"]<-NA
        ##Create name with family and (u.g)
        #res3$gs <- ifelse(is.na(res3$Species),paste0(res3$Genus,"_",res3$ASV), paste0(res3$Genus,"_",res3$Species,"_",res3$ASV))
        #res3$gs <- ifelse(is.na(res3$Genus),paste0(res3$Family, "_",res3$ASV), res3$gs)
        #res3$gs <- ifelse(is.na(res3$Family), paste0(res3$Order,"_",res3$ASV),res3$gs)
        #res3$gs <- ifelse(is.na(res3$Order), paste0(res3$Class, "_",res3$ASV),res3$gs)
        #res3$gs <- ifelse(is.na(res3$Class), paste0(res3$Phylum,"_",res3$ASV),res3$gs)
        ##Make the full trail the First Column
        #res3$names <- res3$ASV
        #res3$Gene.symbol <- res3$gs    
        
        #Drop unnecessary Variables
        drops <- c("Domain","Phylum","Class","Order","Family","Genus","OTU","gs","Species","ASV","Kingdom","names","row2")
        res <- res[ , !(names(res) %in% drops)]
        #res1 <- res1[ , !(names(res1) %in% drops)]
        #res2 <- res2[ , !(names(res2) %in% drops)]
        #res3 <- res3[ , !(names(res3) %in% drops)]
#
        ##merge
        #res5 <- rbind(res1,res2,res3)
        ##only signficant
        #res5 <- res5[!is.na(res5$padj),]
        #res5 <- res5[res5$padj<0.2,]
        #wanted<-which(rownames(res4) %in% rownames(res5))
        ##keep only matching IDs from count data
        #res6<-res4[wanted,]
        #write.table(res6,file="Significant.Taxa.BAL.COPD.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        #meta <- meta %>% select(Study_Linked_ID,three_groups)
        #write.table(meta,file="Significant.Taxa.BAL.COPD.metadata.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

        #----------------------
        ##TABLES
        #----------------------
        #Get Assay Data For Compairson 1
        #GenusData <-as.data.frame(assay(ddsvbal1)) #pruned to selected Genuses based on abundance
        coldata2 <- coldata3
        
        GenusData <-as.data.frame(assay(ddsvbal)) #pruned to selected Genuses based on abundance
        #df <- d3
        #Create Relative Abundance Table
        df <-
            GenusData %>% 
                rownames_to_column('gs') %>%
                group_by(gs) %>% 
                summarise_all(funs(sum)) %>%
                mutate_if(is.numeric, funs(./sum(.))) %>%
                column_to_rownames('gs')
        #Get the ColData for Each Comparison
        coldata.1 <- coldata2[coldata2$Subject_Type=="Smoker.Control",] %>%
                    select(Study_Linked_ID)
        coldata.2 <- coldata2[coldata2$Subject_Type=="COPD",] %>%
                    select(Study_Linked_ID)
        #keep Count data only for each comparison
        needed<-which(colnames(df) %in% rownames(coldata.1))    
        df.1 <- df[,needed]
        needed2<-which(colnames(df) %in% rownames(coldata.2))    
        df.2 <- df[,needed2]
        #Convert Resuts table into a data.frame
        res <- as.data.frame(res)
        #decide what otu to save 
        otu.to.save <-as.character(rownames(res))
        #from relative table we should get the mean across the row of the otu table
        df.1.meanRA <- rowMeans(df.1)
        df.2.meanRA <- rowMeans(df.2)
        #need to subset AND reorder just the otus that we have 
        df.1.meanRA.save <- df.1.meanRA[otu.to.save]
        df.2.meanRA.save <- df.2.meanRA[otu.to.save]
        #add the abundnace data for the res dataframe
        res$abundance.1 <- df.1.meanRA.save
        res$abundance.2 <- df.2.meanRA.save
        #Set Names of Results Table
        res <- setNames(cbind(rownames(res), res, row.names = NULL), c("ASV","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Gene.symbol","abundance.1","abundance.2")) 

        
        #GenusData <-as.data.frame(assay(ddsvbal1)) #pruned to selected Genuses based on abundance
        ##df <- d3
        ##Create Relative Abundance Table
        #df <-
        #    GenusData %>% 
        #        rownames_to_column('gs') %>%
        #        group_by(gs) %>% 
        #        summarise_all(funs(sum)) %>%
        #        mutate_if(is.numeric, funs(./sum(.))) %>%
        #        column_to_rownames('gs')
        ##Get the ColData for Each Comparison
        #coldata.1 <- coldata2[coldata2$Subject_Type_2=="Asymptomatic_SC",] %>%
        #            select(Study_Linked_ID)
        #coldata.2 <- coldata2[coldata2$Subject_Type_2=="COPD",] %>%
        #            select(Study_Linked_ID)
        ##keep Count data only for each comparison
        #needed<-which(colnames(df) %in% rownames(coldata.1))    
        #df.1 <- df[,needed]
        #needed2<-which(colnames(df) %in% rownames(coldata.2))    
        #df.2 <- df[,needed2]
        ##Convert Resuts table into a data.frame
        #res1 <- as.data.frame(res1)
        ##decide what otu to save 
        #otu.to.save <-as.character(rownames(res1))
        ##from relative table we should get the mean across the row of the otu table
        #df.1.meanRA <- rowMeans(df.1)
        #df.2.meanRA <- rowMeans(df.2)
        ##need to subset AND reorder just the otus that we have 
        #df.1.meanRA.save <- df.1.meanRA[otu.to.save]
        #df.2.meanRA.save <- df.2.meanRA[otu.to.save]
        ##add the abundnace data for the res dataframe
        #res1$abundance.1 <- df.1.meanRA.save
        #res1$abundance.2 <- df.2.meanRA.save
        ##Set Names of Results Table
        #res1 <- setNames(cbind(rownames(res1), res1, row.names = NULL), c("ASV","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Gene.symbol","abundance.1","abundance.2")) 
        ##Get Assay Data For Compairson 2
        ##GenusData <-as.data.frame(assay(ddsvbal2)) #pruned to selected Genuses based on abundance
        ##Create Relative Abundance Table
        #GenusData <-as.data.frame(assay(ddsvbal2)) #pruned to selected Genuses based on abundance
        ##Create Relative Abundance Table
        #df <-
        #    GenusData %>% 
        #        rownames_to_column('gs') %>%
        #        group_by(gs) %>% 
        #        summarise_all(funs(sum)) %>%
        #        mutate_if(is.numeric, funs(./sum(.))) %>%
        #        column_to_rownames('gs')
        ##Get the ColData for Each Comparison
        #coldata.1 <- coldata2[coldata2$Subject_Type_2=="Asymptomatic_SC",] %>%
        #            select(Study_Linked_ID)
        #coldata.2 <- coldata2[coldata2$Subject_Type_2=="Symptomatic_SC",] %>%
        #            select(Study_Linked_ID)
        ###keep Count data only for each comparison
        #needed<-which(colnames(df) %in% rownames(coldata.1))    
        #df.1 <- df[,needed]
        #needed2<-which(colnames(df) %in% rownames(coldata.2))    
        #df.2 <- df[,needed2]
        ##Convert Resuts table into a data.frame
        #res2 <- as.data.frame(res2)
        ##decide what otu to save 
        #otu.to.save <-as.character(rownames(res2))
        ##from relative table we should get the mean across the row of the otu table
        #df.1.meanRA <- rowMeans(df.1)
        #df.2.meanRA <- rowMeans(df.2)
        ##need to subset AND reorder just the otus that we have 
        #df.1.meanRA.save <- df.1.meanRA[otu.to.save]
        #df.2.meanRA.save <- df.2.meanRA[otu.to.save]
        ##add the abundnace data for the res dataframe
        #res2$abundance.1 <- df.1.meanRA.save
        #res2$abundance.2 <- df.2.meanRA.save
        ##Set Names of Results Table
        #res2 <- setNames(cbind(rownames(res2), res2, row.names = NULL), c("ASV","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Gene.symbol","abundance.1","abundance.2")) 
        ##Get Assay Data For Compairson 3
        ##GenusData <-as.data.frame(assay(ddsvbal3)) #pruned to selected Genuses based on abundance
        ###Create Relative Abundance Table
        #GenusData <-as.data.frame(assay(ddsvbal3)) #pruned to selected Genuses based on abundance
        ##Create Relative Abundance Table
        #df <-
        #    GenusData %>% 
        #        rownames_to_column('gs') %>%
        #        group_by(gs) %>% 
        #        summarise_all(funs(sum)) %>%
        #        mutate_if(is.numeric, funs(./sum(.))) %>%
        #        column_to_rownames('gs')
        ###Get the ColData for Each Comparison
        #coldata.1 <- coldata2[coldata2$Subject_Type_2=="Symptomatic_SC",] %>%
        #            select(Study_Linked_ID)
        #coldata.2 <- coldata2[coldata2$Subject_Type_2=="COPD",] %>%
        #            select(Study_Linked_ID)
        ###keep Count data only for each comparison
        #needed<-which(colnames(df) %in% rownames(coldata.1))    
        #df.1 <- df[,needed]
        #needed2<-which(colnames(df) %in% rownames(coldata.2))    
        #df.2 <- df[,needed2]
        ##Convert Resuts table into a data.frame
        #res3 <- as.data.frame(res3)
        ##decide what otu to save 
        #otu.to.save <-as.character(rownames(res3))
        ##from relative table we should get the mean across the row of the otu table
        #df.1.meanRA <- rowMeans(df.1)
        #df.2.meanRA <- rowMeans(df.2)
        ##need to subset AND reorder just the otus that we have 
        #df.1.meanRA.save <- df.1.meanRA[otu.to.save]
        #df.2.meanRA.save <- df.2.meanRA[otu.to.save]
        ##add the abundnace data for the res dataframe
        #res3$abundance.1 <- df.1.meanRA.save
        #res3$abundance.2 <- df.2.meanRA.save
        ##Set Names of Results Table
        #res3 <- setNames(cbind(rownames(res3), res3, row.names = NULL), c("ASV","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Gene.symbol","abundance.1","abundance.2")) 
        ##Write Tables of Differential Analysis
        #write.table(res1,file=paste0(counts,".Asymptomatic_SC_vs_COPD_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        #write.table(res2,file=paste0(counts,".Asymptomatic_SC_vs_Sympomatic_SC_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        #write.table(res3,file=paste0(counts,".Sympomatic_SC_vs_COPD_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)        
        
        #Remove Any Data without LOGFC data
        res <- res[!is.na(res$logFC),]
        #res1 <- res1[!is.na(res1$logFC),]
        #res2 <- res2[!is.na(res2$logFC),]
        #res3 <- res3[!is.na(res3$logFC),]
        ##Match the three tables with the same Taxa
        #needed<-which(res1$Gene.symbol %in% res3$Gene.symbol)    
        #res1 <- res1[needed,]
        #needed<-which(res2$Gene.symbol %in% res1$Gene.symbol)  
        #res2 <- res2[needed,]
        #needed<-which(res3$Gene.symbol %in% res1$Gene.symbol)  
        #res3 <- res3[needed,]
        #needed<-which(res1$Gene.symbol %in% res2$Gene.symbol)    
        #res1 <- res1[needed,]
        #needed<-which(res3$Gene.symbol %in% res2$Gene.symbol)    
        #res3 <- res3[needed,]
        # Reorder Results based on FDR for comparison 1
        res = res[order(res$adj.P.Val, na.last = TRUE), ]
        # Keep only top 10
        res <- res %>% slice(1:25)
        #Create Order
        res <- res %>% mutate(start = 1:n())
        # Reorder Results based on FDR for comparison 1
        #res1 = res1[order(res1$adj.P.Val, na.last = TRUE), ]
        ## Keep only top 10
        #res1order <- res1 %>% slice(1:10)
        ## Reorder Results based on FDR for comparison 2
        #res2 = res2[order(res2$adj.P.Val, na.last = TRUE), ]
        ## Keep only top 10
        #res2order <- res2 %>% slice(1:10)
        ## Reorder Results based on FDR for comparison 3
        #res3 = res3[order(res3$adj.P.Val, na.last = TRUE), ]
        ## Keep only top 10
        #res3order <- res3 %>% slice(1:10)
        ##Combine order of three tables
        #resorder <- rbind(res1order,res2order,res3order)
        ##Remove any duplicates
        #resorder <- resorder[!duplicated(resorder$Gene.symbol), ]
        ##Keep only matching Taxa for all three tables
        #res1 <- res1[res1$Gene.symbol %in% resorder$Gene.symbol,]
        ##Set The order
        #res1 <- res1[ order(match(res1$Gene.symbol, resorder$Gene.symbol)), ]
        ##Create Variable for order
        #res1 <- res1 %>% mutate(start = 1:n())
        ##Keep only matching Taxa for all three tables
        #res2 <- res2[res2$Gene.symbol %in% resorder$Gene.symbol,]
        ##Set The order
        #res2 <- res2[ order(match(res2$Gene.symbol, resorder$Gene.symbol)), ]
        ##Keep only matching Taxa for all three tables
        #res2 <- res2 %>% mutate(start = 1:n())
        ##Keep only matching Taxa for all three tables
        #res3 <- res3[res3$Gene.symbol %in% resorder$Gene.symbol,]
        ##Set The order
        #res3 <- res3[ order(match(res3$Gene.symbol, resorder$Gene.symbol)), ]
        ##Keep only matching Taxa for all three tables
        #res3 <- res3 %>% mutate(start = 1:n())
        ##Set Variable for the three comparisson
        #res1$group <- "A"
        #res2$group <- "B"
        #res3$group <- "C"
        #Convert Important columns to Numeric
        res$adj.P.Val <-   as.numeric(as.character(res$adj.P.Val))
        res$logFC <-       as.numeric(as.character(res$logFC))
        res$abundance.1 <- as.numeric(as.character(res$abundance.1))
        res$abundance.2 <- as.numeric(as.character(res$abundance.2))
        #res1$adj.P.Val <-   as.numeric(as.character(res1$adj.P.Val))
        #res1$logFC <-       as.numeric(as.character(res1$logFC))
        #res1$abundance.1 <- as.numeric(as.character(res1$abundance.1))
        #res1$abundance.2 <- as.numeric(as.character(res1$abundance.2))
        #res2$adj.P.Val <-   as.numeric(as.character(res2$adj.P.Val))
        #res2$logFC <-       as.numeric(as.character(res2$logFC))
        #res2$abundance.1 <- as.numeric(as.character(res2$abundance.1))
        #res2$abundance.2 <- as.numeric(as.character(res2$abundance.2))
        #res3$adj.P.Val <-   as.numeric(as.character(res3$adj.P.Val))
        #res3$logFC <-       as.numeric(as.character(res3$logFC))
        #res3$abundance.1 <- as.numeric(as.character(res3$abundance.1))
        #res3$abundance.2 <- as.numeric(as.character(res3$abundance.2))
        ##Bind the three Comparison TAbles
        #resy <- rbind(res1,res2,res3)
        #Replace NA
        res <- res %>% mutate(adj.P.Val = if_else(is.na(adj.P.Val), 0.9, adj.P.Val))
        #resy <- resy %>% mutate(adj.P.Val = if_else(is.na(adj.P.Val), 0.9, adj.P.Val))
        #Create Variable for Color based on Comparison, FDR and LOGFC
        res$col <- ifelse(res$adj.P.Val<0.2 & res$logFC>0, "A",
                   ifelse(res$adj.P.Val<0.2 & res$logFC<0, "B","D"))

        #resy$col <- ifelse(resy$group=="A" & resy$adj.P.Val<0.2 & resy$logFC>0, "A",
        #           ifelse(resy$group=="A" & resy$adj.P.Val<0.2 & resy$logFC<0, "B",
        #           ifelse(resy$group=="B" & resy$adj.P.Val<0.2 & resy$logFC>0, "C",
        #           ifelse(resy$group=="B" & resy$adj.P.Val<0.2 & resy$logFC<0, "B",
        #           ifelse(resy$group=="C" & resy$adj.P.Val<0.2 & resy$logFC>0, "A",
        #           ifelse(resy$group=="C" & resy$adj.P.Val<0.2 & resy$logFC<0, "C","D"))))))
       # 
        #d2
        #results = cbind(as(d2, "data.frame"), as(tax_table(BAL.OTU.Table)[rownames(d2), ], "matrix"))
        #results <- as.data.frame(results)
       # 
        ##convert to character
        #results$ASV <- rownames(results)
        ##Replace any no genus annotation as NA
        #results[results==" s__"]<-NA
        #results[results==" g__"]<-NA
        #results[results==" f__"]<-NA
        #results[results==" o__"]<-NA
        #results[results==" c__"]<-NA
        #results[results=="s__"]<-NA
        #results[results=="g__"]<-NA
        #results[results=="f__"]<-NA
        #results[results=="o__"]<-NA
        #results[results=="c__"]<-NA
        ##Create name with family and (u.g)
        #results$gs <- ifelse(is.na(results$Species),paste0(results$Genus,"_",results$ASV), paste0(results$Genus,"_",results$Species,"_",results$ASV))
        #results$gs <- ifelse(is.na(results$Genus),paste0(results$Family, "_",results$ASV), results$gs)
        #results$gs <- ifelse(is.na(results$Family), paste0(results$Order,"_",results$ASV),results$gs)
        #results$gs <- ifelse(is.na(results$Order), paste0(results$Class, "_",results$ASV),results$gs)
        #results$gs <- ifelse(is.na(results$Class), paste0(results$Phylum,"_",results$ASV),results$gs)
        ##Make the full trail the First Column
        #results$names <- results$ASV
        #results$Gene.symbol <- results$gs
        #results2 <- results
        #resy2 <- resy
        #rownames(results2) <- results2$Gene.symbol
        #resy2 <- resy2[!duplicated(resy2$Gene.symbol), ]

        #rownames(resy2) <- resy2$Gene.symbol
        #wanted<-which(rownames(results2) %in% rownames(resy2))
        ##keep only matching IDs from count data
        #results2<-results2[wanted,]
        #results2 <- results2 %>% select(-c(gs,names,Gene.symbol))
        ##get only BAL
        #balmeta <- colData(ddsvbal)
        #trail <- results2 %>% select(Kingdom,Phylum,Class,Order,Family,Genus,Species,ASV)
        #results3 <- results2 %>% select(-c(Kingdom,Phylum,Class,Order,Family,Genus,Species,ASV))
        #wanted <- which(colnames(results3) %in% rownames(balmeta))
        #results3 <- results3[,wanted]
        #results3 <- cbind(results3,trail)
        #write.table(results3,file="Significant.Taxa.BAL.COPD.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        ##Merge with contam data
        #contam3 <- merge(res1,contam2,by="Gene.symbol",all.x=TRUE)
        #contam3 <- contam3[ order(match(contam3$Gene.symbol, resorder$Gene.symbol)), ]
        ##Create Color vector for tick labels
        #colvec <- ifelse(contam3$contaminant==1, "red", "black")
        #Remove g__
        res$Gene.symbol <- gsub("g__","", res$Gene.symbol,fixed = TRUE)   
        res$Gene.symbol <- gsub("o__","", res$Gene.symbol,fixed = TRUE)   
        res$Gene.symbol <- gsub("f__","", res$Gene.symbol,fixed = TRUE)   
        res$Gene.symbol <- gsub("s__","_",res$Gene.symbol,fixed = TRUE)   


        #resy$Gene.symbol <- gsub("g__","",resy$Gene.symbol,fixed = TRUE)   
        #resy$Gene.symbol <- gsub("o__","",resy$Gene.symbol,fixed = TRUE)   
        #resy$Gene.symbol <- gsub("f__","",resy$Gene.symbol,fixed = TRUE)   
        #resy$Gene.symbol <- gsub("s__","_",resy$Gene.symbol,fixed = TRUE)   
        ##convert abundance to numeric
        res$abundance.2 <- as.numeric(as.character(res$abundance.2))
        res$abundance.1 <- as.numeric(as.character(res$abundance.1))
        res$abundance <- ifelse(res$adj.P.Val<0.2 & res$logFC>0, res$abundance.2, ifelse(res$adj.P.Val<0.2 & res$logFC<0, res$abundance.1,0))
        #resy$abundance.2 <- as.numeric(as.character(resy$abundance.2))
        #resy$abundance.1 <- as.numeric(as.character(resy$abundance.1))
        #resy$abundance <- ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, resy$abundance.1,0))        
        
        #PLOT IT
            ggsave(filename=paste0(counts,".BAL_two_groups_DESEQ2.pdf"),
           ggplot(res, aes(y=reorder(Gene.symbol,-start), x=logFC,fill=col,size=abundance)) +
           #facet_grid(~ group, scales = "free_y")+
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 500 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 500 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Bacterial Size
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 5000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 5000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #DNA Bacterial Size
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 1000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 1000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Virome Size
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 1000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 1000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Fungi Size
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 50000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 50000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+
           geom_point(color="black",alpha=0.8,shape=21)+
           #geom_segment(aes(x=0, y=reorder(Gene.symbol,-start), xend=logFC, yend=reorder(Gene.symbol,-start)), color= "black")+ 
           #geom_segment(data=resy[resy$adj.P.Val<0.2,],aes(yend=reorder(Gene.symbol,-start)), xend=(-30), color= ifelse(resy$adj.P.Val<0.2,"black",NULL), linetype= ifelse(resy$adj.P.Val<0.2,"solid","dashed"))+ 
           geom_segment(data=res[res$adj.P.Val<0.2,],aes(yend=reorder(Gene.symbol,-start)), xend=(-30), color= "black", linetype = "solid",size=1)+ 
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 200000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 200000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Phage size
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 200000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 200000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #DNA Fungi Size
           #scale_fill_manual(values=c("#D01C8B","white"))+ #RNA Virome Colors
           #scale_fill_manual(values=c("#D01C8B","#4DAC26","#FFD479","white"))+ #RNA Fungi Colors/Bacterial Colors
           scale_fill_manual(values=c("B"="#4DAC26","C"="#FFD479","A"="#D01C8B","D"="white"))+ #RNA Phages Colors
           #scale_fill_manual(values=c("#FFD479","white"))+ 
           #scale_fill_manual(values=c("#D01C8B","white"))+ #DNA Fungi Colors
           #scale_fill_manual(values=c("white"))+ #DNA Virome Colors/DNA Phage Colors
           scale_size_continuous(range=c(5, 20),guide=FALSE)+
           theme(panel.background = element_blank(),
               panel.border=element_rect(fill=NA),
               panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
               panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.title=element_text(size=20,face="bold"),
               axis.text.x=element_text(colour="black", size=18, face="bold"),
               axis.text.y=element_text(colour="black",face="bold",size=10),
               axis.ticks=element_line(colour="black"),
        		legend.background = element_rect(color=NA))+
           xlab("") +
           ylab("")+
           #xlim(-7,7)+
           geom_vline(xintercept=0, color="red",linetype="dashed")+
           guides(fill=FALSE),
           width=10, height=10)
           #Create Plot To Get Legened
       p <-ggplot(res, aes(y=reorder(Gene.symbol,-start), x=logFC,fill=col,size=abundance)) +
           #facet_grid(~ group, scales = "free_y")+
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 500 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 500 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Bacterial Size
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 5000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 5000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #DNA Bacterial Size
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 1000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 1000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Virome Size
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 1000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 1000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Fungi Size
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 50000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 50000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+
           geom_point(color="black",alpha=0.8,shape=21)+
           geom_segment(data=res[res$adj.P.Val<0.2,],aes(yend=reorder(Gene.symbol,-start)), xend=(-30), color= "black", linetype = "solid",size=1)+ 
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 200000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 200000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Phage size
           #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 200000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 200000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #DNA Fungi Size
           #scale_fill_manual(values=c("#D01C8B","white"))+ #RNA Virome Colors
           #scale_fill_manual(values=c("#D01C8B","#4DAC26","#FFD479","white"))+ #RNA Fungi Colors/Bacterial Colors
           scale_fill_manual(values=c("B"="#4DAC26","C"="#FFD479","A"="#D01C8B","D"="white"))+ #RNA Phages Colors
           #scale_fill_manual(values=c("#FFD479","white"))+ 
           #scale_fill_manual(values=c("#D01C8B","white"))+ #DNA Fungi Colors
           #scale_fill_manual(values=c("white"))+ #DNA Virome Colors/DNA Phage Colors
           #scale_size_area(max_size=20)+
           #scale_radius(range=c(1, 6))+
           scale_size_continuous(name="Relative Abundance",range=c(5, 20))+
           theme(panel.background = element_blank(),
               panel.border=element_rect(fill=NA),
               panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
               panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.title=element_text(size=20,face="bold"),
               axis.text.x=element_text(colour="black", size=18, face="bold"),
               axis.text.y=element_text(colour="black",face="bold",size=10),
               axis.ticks=element_line(colour="black"),
        	   legend.background = element_rect(color=NA),
               legend.key = element_rect(colour = "transparent", fill = "white"))+
           xlab("") +
           ylab("")+
           #xlim(-7,7)+
           geom_vline(xintercept=0, color="red",linetype="dashed")+
           guides(fill=FALSE)  
        # Extract the legend. Returns a gtable
            leg <- get_legend(p)
        # Convert to a ggplot and print
           ggsave(filename=paste0(counts,".BAL_two_groups_DESEQ2_legend.pdf"),
            as_ggplot(leg),
           width=5, height=5)
        
        ##PLOT IT
        #   ggsave(filename=paste0(counts,".BAL_three_groups_DESEQ2.pdf"),
        #   ggplot(resy, aes(y=reorder(Gene.symbol,-start), x=logFC,fill=col,size=abundance)) +
        #   facet_grid(~ group, scales = "free_y")+
        #   #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 500 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 500 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Bacterial Size
        #   #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 5000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 5000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #DNA Bacterial Size
        #   #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 1000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 1000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Virome Size
        #   #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 1000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 1000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Fungi Size
        #   # geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 50000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 50000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+
        #   geom_point(color="black",alpha=0.8,shape=21)+
        #   geom_segment(data=resy[resy$adj.P.Val<0.2,],aes(yend=reorder(Gene.symbol,-start)), xend=(-30), color= "black", linetype = "solid",size=1)+ 
        #   #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 200000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 200000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Phage size
        #   #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 200000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 200000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #DNA Fungi Size
        #   #scale_fill_manual(values=c("#D01C8B","white"))+ #RNA Virome Colors
        #   #scale_fill_manual(values=c("#D01C8B","#4DAC26","#FFD479","white"))+ #RNA Fungi Colors/Bacterial Colors
        #   scale_fill_manual(values=c("B"="#4DAC26","C"="#FFD479","A"="#D01C8B","D"="white"))+ #RNA Phages Colors
        #   #scale_fill_manual(values=c("#FFD479","white"))+ 
        #   #scale_fill_manual(values=c("#D01C8B","white"))+ #DNA Fungi Colors
        #   #scale_fill_manual(values=c("white"))+ #DNA Virome Colors/DNA Phage Colors
        #   scale_size_continuous(range=c(5, 20),guide=FALSE)+
        #   theme(panel.background = element_blank(),
        #       panel.border=element_rect(fill=NA),
        #       panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        #       panel.grid.minor = element_blank(),
        #       strip.background=element_blank(),
        #       axis.title=element_text(size=20,face="bold"),
        #       axis.text.x=element_text(colour="black", size=18, face="bold"),
        #       axis.text.y=element_text(colour="black",face="bold",size=10),
        #       axis.ticks=element_line(colour="black"),
        #		legend.background = element_rect(color=NA))+
        #   xlab("") +
        #   ylab("")+
        #   #xlim(-7,7)+
        #	geom_vline(xintercept=0, color="red",linetype="dashed")+
        #	guides(fill=FALSE),
        #   width=20, height=5)
        #   #Create Plot To Get Legened
       #p <-ggplot(resy, aes(y=reorder(Gene.symbol,-start), x=logFC,fill=col,size=abundance)) +
       #    facet_grid(~ group, scales = "free_y")+
       #    #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 500 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 500 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Bacterial Size
       #    #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 5000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 5000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #DNA Bacterial Size
       #    #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 1000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 1000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Virome Size
       #    #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 1000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 1000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Fungi Size
       #    #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 50000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 50000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+
       #    #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+
       #    geom_point(color="black",alpha=0.8,shape=21)+
       #    geom_segment(data=resy[resy$adj.P.Val<0.2,],aes(yend=reorder(Gene.symbol,-start)), xend=(-30), color= "black", linetype = "solid",size=1)+ 
       #    #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 200000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 200000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Phage size
       #    #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 200000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 200000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #DNA Fungi Size
       #    #scale_fill_manual(values=c("#D01C8B","white"))+ #RNA Virome Colors
       #    #scale_fill_manual(values=c("#D01C8B","#4DAC26","#FFD479","white"))+ #RNA Fungi Colors/Bacterial Colors
       #    scale_fill_manual(values=c("B"="#4DAC26","C"="#FFD479","A"="#D01C8B","D"="white"))+ #RNA Phages Colors
       #    #scale_fill_manual(values=c("#FFD479","white"))+ 
       #    #scale_fill_manual(values=c("#D01C8B","white"))+ #DNA Fungi Colors
       #    #scale_fill_manual(values=c("white"))+ #DNA Virome Colors/DNA Phage Colors
       #    #scale_size_area(max_size=20)+
       #    #scale_radius(range=c(1, 6))+
       #    scale_size_continuous(name="Relative Abundance",range=c(5, 20))+
       #    theme(panel.background = element_blank(),
       #        panel.border=element_rect(fill=NA),
       #        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
       #        panel.grid.minor = element_blank(),
       #        strip.background=element_blank(),
       #        axis.title=element_text(size=20,face="bold"),
       #        axis.text.x=element_text(colour="black", size=18, face="bold"),
       #        axis.text.y=element_text(colour="black",face="bold",size=10),
       #        axis.ticks=element_line(colour="black"),
       # 	   legend.background = element_rect(color=NA),
       #        legend.key = element_rect(colour = "transparent", fill = "white"))+
       #    xlab("") +
       #    ylab("")+
       #    #xlim(-7,7)+
       #    geom_vline(xintercept=0, color="red",linetype="dashed")+
       #    guides(fill=FALSE)  
       # # Extract the legend. Returns a gtable
       #     leg <- get_legend(p)
       # # Convert to a ggplot and print
       #    ggsave(filename=paste0(counts,".BAL_three_groups_DESEQ2_legend.pdf"),
       #     as_ggplot(leg),
       #    width=5, height=5)

}
