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
#library(egg)


#Set Theme for Figures
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")

#Choose Alpha/FDR
alpha = 0.2

#Create Function For Analysis
deseq <- function(metadata,counts,recounts,sampletype,comparison) {
        #Load MetaData
        coldata <- read.delim2(metadata, sep="\t")
        coldata$Study_Linked_ID <- coldata$X #RNA
        #coldata$Study_Linked_ID <- coldata$GTC_Updated_Name_MetatranscriptomeMH2021 #Metatranscriptome
        #coldata$Study_Linked_ID <- gsub("-","_",coldata$Study_Linked_ID) #Metatranscriptome
        #coldata$Study_Linked_ID <- gsub(".Sup","_Sup",coldata$Study_Linked_ID) #Metatranscriptome
        #coldata$Study_Linked_ID <- gsub("-",".",coldata$Study_Linked_ID) #Metagenome
        coldata <- coldata[!duplicated(coldata$Study_Linked_ID), ]
        coldata <- coldata[coldata$Study_Linked_ID!="n.a", ]
        rownames(coldata) <- coldata$Study_Linked_ID
        coldata$Sample.Type <- coldata$Sample_Description_s
        coldata$Subject_Type <- coldata$GOLD_Post_Final_lns
        #coldata$Sample.Type <- ifelse(coldata$Sample.Type=="BKG","UA",
        #                ifelse(coldata$Sample.Type=="UA","BKG",as.character(coldata$Sample.Type))) #Metagenome
        #Load Count Data
        mycounts  <-read.delim2(paste0(counts,".txt"), sep="\t", row.names=1)
        #colnames(mycounts) <- gsub("\\.","_",colnames(mycounts))
        #mycounts   <- mycounts %>% select(!category) #Metagenome       
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
        ddsv <- DESeqDataSetFromMatrix(countData = d1,
                              colData = coldata2,
                              design= ~ Subject_Type)
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
        #Transforming data - Option #2 is variance stabilizing transformation
        vsdv <- varianceStabilizingTransformation(ddsv)
        #vsdv <- rlog(ddsv, fitType="local",blind=FALSE)
        #----------------------
        ##HEATMAP
        #----------------------
        #filter out genes where there are less than 3 samples with normalized counts greater than or equal to 100.
        idx <- rowSums(relcounts3 >= 0.0012 ) >= 18 #Metatranscriptome
        df <- relcounts3[idx,]
        #cluster Genuses(row)
        GenusData.Bray.dist <-vegdist(df, method = "bray")
        #cluster samples(Col)
        Samples.Bray.dist = vegdist(t(df), method="bray")
        #Set Color Scale for Heatmap
        mypalette <- colorRampPalette(c('#ffffff','#4169E1','#0000CD'))(100)
        #Set Colors for each sample type for HeatMap
        annon_colors= list(Subject_Type=c(Smoker.Control="#D01C8B", COPD="#4DAC26"))
        #Choose lables for Samples
        df2 <- data.frame(Subject_Type = colData(ddsv)[,c("Subject_Type")], row.names = rownames(colData(ddsv)))
        #Plot Heatmap
        ggsave(filename=paste0(counts,".Heatmap_Pruned.pdf"),
            pheatmap(df, cluster_rows=TRUE, show_rownames=TRUE, 
            cluster_cols=TRUE,annotation_col=df2,
            scale="row",
            clustering_distance_rows = GenusData.Bray.dist,
            clustering_distance_cols = Samples.Bray.dist,
            clustering_method="mcquitty",
            gaps_col=50,
            border_color="black",
            color = colorRampPalette(c('#4169E1','#ffffff','#0000CD'))(100),
            annotation_colors=annon_colors[1],legend=FALSE),
            height = 20, width = 20)
        #----------------------
        ##Create graph of Reads
        #----------------------
        #Sum number or reads per sample
        summary <- as.data.frame(rowSums(t(assay(ddsv))))
        #rownames(summary) <- colnames(ddsv)
        #Merge Reads Data with MetaData
        require(data.table)
        Reads <- as.data.frame(merge(x = summary, y = colData(ddsv), by = "row.names", all.x = TRUE))
        #Rename Column of Reads
        colnames(Reads)[colnames(Reads)=="rowSums.t.assay.ddsv..."] <- "Reads"
        #Set Order Of Figure
        Reads$or <-ifelse(Reads$Subject_Type=="Smoker.Control", 1, NA)
        Reads$or <-ifelse(Reads$Subject_Type=="COPD", 2, Reads$or)

        #Create Figure
            ggsave(filename=paste0(counts,".Reads.pdf"),
            ggplot(Reads, aes(x= reorder(Subject_Type, +or), y=Reads, fill=Subject_Type)) + 
            stat_boxplot(geom ='errorbar', width=0.1)+
            geom_boxplot(outlier.shape = NA, width=0.5)+
            geom_jitter(shape=1, position=position_jitter(0.2))+
            scale_fill_manual(values=c("#D01C8B","#4DAC26" )) + 
            scale_x_discrete(labels = c('Smoker.Control','COPD'))+ 
            scale_y_continuous(name="Reads",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
            xlab("Subject_Type")+
            theme,
            height = 7, width = 5)
        #----------------------
        ##DESEQ
        #----------------------
        #Covert Variable to Factor
        coldata2$Subject_Type <- as.factor(coldata2$Subject_Type)
        #Create Deseq object for BAL analysis
        ddsvbal <- DESeqDataSetFromMatrix(countData = d1,
                              colData = coldata2,
                              design= ~ Subject_Type)
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsvbal), 1, gm_mean)        
        #Estimate Factors of DESeq Object
        ddsvbal <- estimateSizeFactors(ddsvbal, geoMeans = geoMeans)
        #Calculate geometric means prior to estimate size factor
        #gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        #geoMeans = apply(counts(ddsvbal1), 1, gm_mean)
        ##Estimate Factors of DESeq Object
        #ddsvbal1 <- estimateSizeFactors(ddsvbal1, geoMeans = geoMeans)
        ##Calculate geometric means prior to estimate size factor
        #gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        #geoMeans = apply(counts(ddsvbal2), 1, gm_mean)
        ##Estimate Factors of DESeq Object
        #ddsvbal2 <- estimateSizeFactors(ddsvbal2, geoMeans = geoMeans)
        ##Calculate geometric means prior to estimate size factor
        #gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        #geoMeans = apply(counts(ddsvbal3), 1, gm_mean)
        ##Estimate Factors of DESeq Object
        #ddsvbal3 <- estimateSizeFactors(ddsvbal3, geoMeans = geoMeans)
        #Variance Stablize the data
        vsdvbal <- varianceStabilizingTransformation(ddsvbal)
        #vsdvbal <- rlog(ddsvbal, fitType="local",blind=FALSE)
        #----------------------
        ##PCOA
        #----------------------
        #Create Distance Matrix
        vsdvbal0 <- ifelse(assay(vsdvbal)<0,0,assay(vsdvbal))
        vegdist   = vegdist(t(vsdvbal0), method="bray")
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
        centroids <- aggregate(cbind(PC1,PC2)~ Subject_Type,data= newResults, mean)
        #Merge the Centroid Data into the PCOA Data
        newResults <- merge(newResults,centroids,by="Subject_Type",suffixes=c("",".centroid"))
        #Create Table for Statistics    
        data.adonis <- data.frame(colData(vsdvbal))
        #Run the Statistics
        samplepermanova <- adonis(vegdist ~ Subject_Type, data.adonis)
        samplepermanova <- as.data.frame(samplepermanova$aov.tab)
        samplepermanova <- samplepermanova$'Pr(>F)'[1]
        #PLOT IT
            ggsave(filename=paste0(counts,".BRAY_vsd_Subject_Type_PERMANOVA_",samplepermanova,".pdf"),
            ggplot(newResults, aes(PC1, PC2,color=Subject_Type)) + # Graph PC1 and PC2
            xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
            ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
            scale_color_manual(values=c("#D01C8B","#4DAC26" )) + 
            geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Subject_Type))+ 
            geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=Subject_Type), size=10) + #BKG DFW Lung.Tissue.In Lung.Tissue.UnIn Mock
            geom_point(data=newResults,aes(color=Subject_Type),size=5,alpha=0.5) + # Set the size of the points
            theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
            panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
            panel.grid.minor = element_blank(),strip.background=element_blank(),
            axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
            axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),
            plot.margin=unit(c(1,1,1,1),"line"), legend.position="none"),
            height = 10, width = 10)
        #----------------------
        ##Alpha Diversity
        #----------------------
        #Calcultes Shannon Diversity
        ddsvbal$Shannon = diversity(vsdvbal0, index = "shannon", MARGIN = 2, base = exp(1))
        #ddsvbal$Shannon = diversity(d3, index = "shannon", MARGIN = 2, base = exp(1))
        #Convert to data frame for ggplot
        shannon = as.data.frame(colData(ddsvbal))
        #Set Order Of Figure
        shannon$or <-ifelse(shannon$Subject_Type=="Smoker.Control", 1,NA)
        shannon$or <-ifelse(shannon$Subject_Type=="COPD",2 ,shannon$or)
        #Make Sure Shannon is Numeric
        shannon$Shannon <- as.numeric(as.character(shannon$Shannon))
        #Make sure Group is factor
        shannon$Subject_Type <- as.factor(shannon$Subject_Type)
        #Check Statistics
        sampleshannon <- kruskal.test(Shannon ~ Subject_Type, data = shannon)
        sampleshannon <- sampleshannon$p.value
        #PLOT IT
            ggsave(filename=paste0(counts,".SHANNON_vsd_Subject_Type_Pvalue_",sampleshannon,".pdf"),
            ggplot(shannon, aes(x= reorder(Subject_Type, +or), y=Shannon, fill=Subject_Type)) + 
            stat_boxplot(geom ='errorbar', width=0.1)+
            geom_boxplot(outlier.shape = NA, width=0.5)+
            geom_jitter(shape=1, position=position_jitter(0.2))+
            scale_fill_manual(values=c("#D01C8B","#4DAC26" )) + 
            #scale_x_discrete(labels = c('<28 Days','>28 Days','Deceased'))+ 
            ylab("Shannon Diversity") + 
            xlab("Subject Type")+
            #scale_y_continuous(breaks = seq(2, 5, .5), limits = c(2, 5)) + #For Fungi
            #scale_y_continuous(breaks = seq(3, 5, .5), limits = c(3, 5)) + #For Fungi
            theme,
            useDingbats=FALSE,
            height = 7, width = 5)
        #----------------------
        ##Differential Analysis
        #----------------------
        #DropLevels
        ddsvbal$Subject_Type <- droplevels(ddsvbal$Subject_Type)
        #ddsvbal1$three_groups <- droplevels(ddsvbal1$three_groups)
        #ddsvbal2$three_groups <- droplevels(ddsvbal2$three_groups)
        #ddsvbal3$three_groups <- droplevels(ddsvbal3$three_groups)
        #Set Reference
        ddsvbal$Subject_Type <- relevel(ddsvbal$Subject_Type, ref ="Smoker.Control")
        #ddsvbal1$three_groups <- relevel(ddsvbal1$three_groups, ref ="Less_Than_28_days_on_vent")
        #ddsvbal2$three_groups <- relevel(ddsvbal2$three_groups, ref ="Less_Than_28_days_on_vent")
        #ddsvbal3$three_groups <- relevel(ddsvbal3$three_groups, ref ="Greater_Than_28_days_on_vent")
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
        ##TABLES
        #----------------------
        #Get Assay Data For Compairson 1
        #GenusData <-as.data.frame(assay(ddsvbal1)) #pruned to selected Genuses based on abundance
        df <- d2
        #Create Relative Abundance Table
        #df <-
        #    GenusData %>% 
        #        rownames_to_column('gs') %>%
        #        group_by(gs) %>% 
        #        summarise_all(funs(sum)) %>%
        #        mutate_if(is.numeric, funs(./sum(.))) %>%
        #        column_to_rownames('gs')
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
        res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
        #Get Assay Data For Compairson 2
        #GenusData <-as.data.frame(assay(ddsvbal2)) #pruned to selected Genuses based on abundance
        #Create Relative Abundance Table
        #df <-
        #    GenusData %>% 
        #        rownames_to_column('gs') %>%
        #        group_by(gs) %>% 
        #        summarise_all(funs(sum)) %>%
        #        mutate_if(is.numeric, funs(./sum(.))) %>%
        #        column_to_rownames('gs')
        #Get the ColData for Each Comparison
        #coldata.1 <- coldata2[coldata2$three_groups=="Less_Than_28_days_on_vent",] %>%
        #            select(Study_Linked_ID)
        #coldata.2 <- coldata2[coldata2$three_groups=="Greater_Than_28_days_on_vent",] %>%
        #            select(Study_Linked_ID)
        ##keep Count data only for each comparison
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
        #res2 <- setNames(cbind(rownames(res2), res2, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
        ##Get Assay Data For Compairson 3
        ##GenusData <-as.data.frame(assay(ddsvbal3)) #pruned to selected Genuses based on abundance
        ###Create Relative Abundance Table
        ##df <-
        ##    GenusData %>% 
        ##        rownames_to_column('gs') %>%
        ##        group_by(gs) %>% 
        ##        summarise_all(funs(sum)) %>%
        ##        mutate_if(is.numeric, funs(./sum(.))) %>%
        ##        column_to_rownames('gs')
        ##Get the ColData for Each Comparison
        #coldata.1 <- coldata2[coldata2$three_groups=="Greater_Than_28_days_on_vent",] %>%
        #            select(Study_Linked_ID)
        #coldata.2 <- coldata2[coldata2$three_groups=="Dead",] %>%
        #            select(Study_Linked_ID)
        ##keep Count data only for each comparison
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
        #res3 <- setNames(cbind(rownames(res3), res3, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
        #Write Tables of Differential Analysis
        write.table(res,file=paste0(counts,".Smoker_vs_COPD_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        #write.table(res1,file=paste0(counts,".Less_Than_28D_vs_Dead_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        #write.table(res2,file=paste0(counts,".Less_Than_28D_vs_Greater_Than_28D_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        #write.table(res3,file=paste0(counts,".Greater_Than_28D_vs_Dead_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        res$sig <- -log10(res$adj.P.Val)        
        res[is.infinite(res$sig),"sig"] <- 350
        ## Volcano plot of adjusted p-values
        cols <- densCols(res$logFC, res$sig)
        cols[res$pvalue ==0] <- "purple"
        cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
        cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "green"
        res$pch <- 19
        res$pch[res$pvalue ==0] <- 6
        ggsave(filename=paste0(counts,".BAL_Smoker_vs_COPD_DESEQ2.pdf"),
            ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
            geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 5000 * res$abundance.2, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 5000 * res$abundance.1,1)),alpha=0.6) + #Chose Colors for dots
            geom_text_repel(aes(label=ifelse(res$logFC<(-1) & res$adj.P.Val < alpha , as.character(res$Gene.symbol),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
            geom_text_repel(aes(label=ifelse(res$logFC>1 & res$adj.P.Val < alpha , as.character(res$Gene.symbol),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
            theme(legend.position = "none") +
            geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
            xlab("Effect size: log2(fold-change)") +
            ylab("-log10(adjusted p-value)") + 
            #ylim(0,20)+
            theme,
            width=5, height=5)
        #CLUSTERING
        idx <- rowSums(relcounts3 >= 0.0001 ) >= 55 #Metagenome All
        df <- relcounts3[idx,] #Metagenome All
        #cluster Genuses(row)
        umap_model_1 = umap(as.matrix(t(df)))
        umapcord <- as.data.frame(umap_model_1$layout)
        umapcord %<>% mutate(PC1 = umap_model_1$layout[, 1], PC2 = umap_model_1$layout[,2])
        rownames(umapcord) <- rownames(coldata2)
        tsne_model_1 = Rtsne(as.matrix(t(df)), check_duplicates=FALSE, pca=TRUE,theta=0.5,perplexity=5, dims=2)
        tsnecord <- as.data.frame(tsne_model_1$Y)
        tsnecord %<>% mutate(PC1 = tsne_model_1$Y[, 1], PC2 = tsne_model_1$Y[,2])
        rownames(tsnecord) <- rownames(coldata2)
        #TSNE MODEL 2
        tsne_model_2 = Rtsne(as.matrix(df), check_duplicates=FALSE, pca=TRUE,theta=0.5,perplexity=5, dims=2)
        tsnecord2 <- as.data.frame(tsne_model_2$Y)
        tsnecord2 %<>% mutate(PC1 = tsne_model_2$Y[, 1], PC2 = tsne_model_2$Y[,2])
        rownames(tsnecord2) <- rownames(df)
        tsnecord2$Sample.Type <- "Taxa"
        tsnecord2$Strain <- rownames(tsnecord2)
        #Clustering
        km.norm.tsne = kmeans(tsne_model_1$Y, 3, nstart = 25)
        km.norm.umap = kmeans(umap_model_1$layout, 3, nstart = 25)
        hc.norm.tsne = hclust(dist(tsne_model_1$Y)) #Complete
        hc.norm.umap = hclust(dist(umap_model_1$layout)) #Complete
        #Put the cluster classification in the PC table
        tsnecord$kmeans = factor(km.norm.tsne$cluster)
        umapcord$kmeans = factor(km.norm.umap$cluster)
        tsnecord$hclust = factor(cutree(hc.norm.tsne, 3))
        umapcord$hclust = factor(cutree(hc.norm.umap, 3))
        #Merge TSNE with MetaData
            tsne <- merge(x = tsnecord, y = coldata2, by = "row.names", all.x = TRUE)
            tsne$col <- ifelse(tsne$kmeans==1,"Cluster_3", 
                       ifelse(tsne$kmeans==2,"Cluster_2","Cluster_1"))
            tsne$col2 <- ifelse(tsne$hclust==1,"Cluster_3", 
                       ifelse(tsne$hclust==2,"Cluster_2","Cluster_1"))            
            #Calculate the Centroid Value
            centroids <- aggregate(cbind(PC1,PC2)~ Sample.Type,data= tsne, mean)
            #Merge the Centroid Data into the PCOA Data
            tsne <- merge(tsne,centroids,by="Sample.Type",suffixes=c("",".centroid"))
            ggsave(filename=paste0(counts,".TSNE_Kmeans.PCOA.pdf"),
                ggplot(tsne, aes(PC1, PC2,color=Sample.Type)) + # Graph PC1 and PC2
                scale_color_manual(values=c("#1700F5", "#BEBEBE", "#F7A501")) + 
                geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Sample.Type))+ 
                geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=Sample.Type), size=10) + #BKG DFW Lung.Tissue.In Lung.Tissue.UnIn Mock
                geom_point(data=tsne,aes(fill=col),color="black",size=6, pch=21,alpha=0.5) + # Set the size of the points
                scale_fill_manual(values=c("#4DAC26","#FFD479","#D01C8B")) + 
                xlab("tsne 1")+
                ylab("tsne 2")+
                theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
                panel.grid.minor = element_blank(),strip.background=element_blank(),
                axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
                axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),
                plot.margin=unit(c(1,1,1,1),"line"),legend.position="none"),
            height = 10, width = 10)
        write.table(tsne,file=paste0(counts,".TSNE_Clustering.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE) 
        #HEATMAP
        idx <- rowSums(relcounts3 >= 0.0012 ) >= 18 #Metagenome All
        heatmap <- relcounts3[idx,] #Metagenome All
        #cluster Genuses(row)
        GenusData.Bray.dist <-vegdist(heatmap, method = "bray")
        #cluster samples(Col)
        Samples.Bray.dist = vegdist(t(heatmap), method="bray")
        #Set Color Scale for Heatmap
        mypalette <- colorRampPalette(c('#ffffff','#4169E1','#0000CD'))(100)
        #Set Colors for each sample type for HeatMap
        annon_colors= list(Subject_Type=c(Smoker.Control="#D01C8B", COPD="#4DAC26"))
        #Choose lables for Samples
        df2 <- data.frame(Subject_Type = tsne$Subject_Type, Cluster = tsne$col, row.names = tsne$Row.names)
        #Plot Heatmap
        ggsave(filename=paste0(counts,".Heatmap_Pruned_with_Clusters.pdf"),
            pheatmap(heatmap, cluster_rows=TRUE, show_rownames=TRUE, 
            cluster_cols=TRUE,annotation_col=df2,
            scale="row",
            clustering_distance_rows = GenusData.Bray.dist,
            clustering_distance_cols = Samples.Bray.dist,
            clustering_method="mcquitty",
            gaps_col=50,
            border_color="black",
            color = colorRampPalette(c('#4169E1','#ffffff','#0000CD'))(100),
            annotation_colors=annon_colors,legend=FALSE),
            height = 20, width = 20)
        #Subset Rel Table for BAL only
        wanted<-which(rownames(d1) %in% rownames(df))
        #keep only matching IDs from count data
        d3<-d1[wanted,]
        #Covert Variable to Factor
        tsne$three_groups <- as.factor(tsne$col)
        #Create Deseq object for BAL analysis
        ddsvbal <- DESeqDataSetFromMatrix(countData = d3,
                              colData = tsne,
                              design= ~ three_groups)
        #Create Seperate tables with each of the different compairsons
        ddsvbal1 <- ddsvbal[, ddsvbal$three_groups %in% c("Cluster_1","Cluster_3")]
        ddsvbal2 <- ddsvbal[, ddsvbal$three_groups %in% c("Cluster_1","Cluster_2")]
        ddsvbal3 <- ddsvbal[, ddsvbal$three_groups %in% c("Cluster_2","Cluster_3")]
        #Create Tables of First Comparison
        balcounts1 <- assay(ddsvbal1)
        balmeta1   <- colData(ddsvbal1)
        #Create Tables of Second Comparison
        balcounts2 <- assay(ddsvbal2)
        balmeta2   <- colData(ddsvbal2)
        #Create Tables of Third Comparison
        balcounts3 <- assay(ddsvbal3)
        balmeta3   <- colData(ddsvbal3)
        #Create new DESEQ2 Objects
        ddsvbal1 <- DESeqDataSetFromMatrix(countData = balcounts1,
                                      colData = balmeta1,
                                      design= ~ three_groups)
        ddsvbal2 <- DESeqDataSetFromMatrix(countData = balcounts2,
                                      colData = balmeta2,
                                      design= ~ three_groups)
        ddsvbal3 <- DESeqDataSetFromMatrix(countData = balcounts3,
                                    colData = balmeta3,
                                    design= ~ three_groups)
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsvbal), 1, gm_mean)        
        #Estimate Factors of DESeq Object
        ddsvbal <- estimateSizeFactors(ddsvbal, geoMeans = geoMeans)
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsvbal1), 1, gm_mean)
        #Estimate Factors of DESeq Object
        ddsvbal1 <- estimateSizeFactors(ddsvbal1, geoMeans = geoMeans)
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsvbal2), 1, gm_mean)
        #Estimate Factors of DESeq Object
        ddsvbal2 <- estimateSizeFactors(ddsvbal2, geoMeans = geoMeans)
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsvbal3), 1, gm_mean)
        #Estimate Factors of DESeq Object
        ddsvbal3 <- estimateSizeFactors(ddsvbal3, geoMeans = geoMeans)
        #Variance Stablize the data
        #vsdvbal <- varianceStabilizingTransformation(ddsvbal)
        #vsdvbal <- rlog(ddsvbal, fitType="local",blind=FALSE)
        #----------------------
        ##Alpha Diversity
        #----------------------
        #Calcultes Shannon Diversity
        ddsvbal$Shannon = diversity(assay(ddsvbal), index = "shannon", MARGIN = 2, base = exp(1))
        #ddsvbal$Shannon = diversity(d3, index = "shannon", MARGIN = 2, base = exp(1))
        #Convert to data frame for ggplot
        shannon = as.data.frame(colData(ddsvbal))
        #Set Order Of Figure
        shannon$or <-ifelse(shannon$three_groups=="Cluster_1", 1,NA)
        shannon$or <-ifelse(shannon$three_groups=="Cluster_2",2 ,shannon$or)
        shannon$or <-ifelse(shannon$three_groups=="Cluster_3",3,shannon$or)
        #Make Sure Shannon is Numeric
        shannon$Shannon <- as.numeric(as.character(shannon$Shannon))
        #Make sure Group is factor
        shannon$three_groups <- as.factor(shannon$three_groups)
        #Check Statistics
        sampleshannon <- kruskal.test(Shannon ~ three_groups, data = shannon)
        sampleshannon <- sampleshannon$p.value
        #PLOT IT
            ggsave(filename=paste0(counts,".SHANNON_vsd_Clusters_Pvalue_",sampleshannon,".pdf"),
            ggplot(shannon, aes(x= reorder(three_groups, +or), y=Shannon, fill=three_groups)) + 
            stat_boxplot(geom ='errorbar', width=0.1)+
            geom_boxplot(outlier.shape = NA, width=0.5)+
            geom_jitter(shape=1, position=position_jitter(0.2))+
            scale_fill_manual(values=c("#4DAC26","#FFD479","#D01C8B")) + 
            scale_x_discrete(labels = c('Cluster 1','Cluster 2','Cluster 3'))+ 
            ylab("Shannon Diversity") + 
            xlab("Clusters")+
            #scale_y_continuous(breaks = seq(2, 5, .5), limits = c(2, 5)) + #For Fungi
            scale_y_continuous(breaks = seq(3, 5, .5), limits = c(3, 5)) + #For Fungi
            theme,
            useDingbats=FALSE,
            height = 7, width = 5)
        #Check Statistics
        kruskal.test(Shannon ~ three_groups, data = shannon)
        #----------------------
        ##Differential Analysis
        #----------------------
        #DropLevels
        ddsvbal1$three_groups <- droplevels(ddsvbal1$three_groups)
        ddsvbal2$three_groups <- droplevels(ddsvbal2$three_groups)
        ddsvbal3$three_groups <- droplevels(ddsvbal3$three_groups)
        #Set Reference
        ddsvbal1$three_groups <- relevel(ddsvbal1$three_groups, ref ="Cluster_1")
        ddsvbal2$three_groups <- relevel(ddsvbal2$three_groups, ref ="Cluster_1")
        ddsvbal3$three_groups <- relevel(ddsvbal3$three_groups, ref ="Cluster_2")
        #Run the differential Analysis: Lung Cancer Vs Wild Type --> positive is upregulated in Lung Cancer; Negative is down regulated
        ddsvbal1  <- DESeq(ddsvbal1)
        ddsvbal2  <- DESeq(ddsvbal2)
        ddsvbal3  <- DESeq(ddsvbal3)
        #Output Result Table
        res1     <- results(ddsvbal1, cooksCutoff=FALSE)
        res2     <- results(ddsvbal2, cooksCutoff=FALSE)
        res3     <- results(ddsvbal3, cooksCutoff=FALSE)
        #----------------------
        ##TABLES
        #----------------------
        #Get Assay Data For Compairson 1
        #GenusData <-as.data.frame(assay(ddsvbal1)) #pruned to selected Genuses based on abundance
        coldata2 <- as.data.frame(colData(ddsvbal))
        #Create Relative Abundance Table
        #df <-
        #    GenusData %>% 
        #        rownames_to_column('gs') %>%
        #        group_by(gs) %>% 
        #        summarise_all(funs(sum)) %>%
        #        mutate_if(is.numeric, funs(./sum(.))) %>%
        #        column_to_rownames('gs')
        #Get the ColData for Each Comparison
        coldata.1 <- coldata2[coldata2$three_groups=="Cluster_1",] %>%
                    select(Study_Linked_ID)
        coldata.2 <- coldata2[coldata2$three_groups=="Cluster_3",] %>%
                    select(Study_Linked_ID)
        #keep Count data only for each comparison
        needed<-which(colnames(df) %in% rownames(coldata.1))    
        df.1 <- df[,needed]
        needed2<-which(colnames(df) %in% rownames(coldata.2))    
        df.2 <- df[,needed2]
        #Convert Resuts table into a data.frame
        res1 <- as.data.frame(res1)
        #decide what otu to save 
        otu.to.save <-as.character(rownames(res1))
        #from relative table we should get the mean across the row of the otu table
        df.1.meanRA <- rowMeans(df.1)
        df.2.meanRA <- rowMeans(df.2)
        #need to subset AND reorder just the otus that we have 
        df.1.meanRA.save <- df.1.meanRA[otu.to.save]
        df.2.meanRA.save <- df.2.meanRA[otu.to.save]
        #add the abundnace data for the res dataframe
        res1$abundance.1 <- df.1.meanRA.save
        res1$abundance.2 <- df.2.meanRA.save
        #Set Names of Results Table
        res1 <- setNames(cbind(rownames(res1), res1, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
        #Get Assay Data For Compairson 2
        #GenusData <-as.data.frame(assay(ddsvbal2)) #pruned to selected Genuses based on abundance
        #Create Relative Abundance Table
        #df <-
        #    GenusData %>% 
        #        rownames_to_column('gs') %>%
        #        group_by(gs) %>% 
        #        summarise_all(funs(sum)) %>%
        #        mutate_if(is.numeric, funs(./sum(.))) %>%
        #        column_to_rownames('gs')
        #Get the ColData for Each Comparison
        coldata.1 <- coldata2[coldata2$three_groups=="Cluster_1",] %>%
                    select(Study_Linked_ID)
        coldata.2 <- coldata2[coldata2$three_groups=="Cluster_2",] %>%
                    select(Study_Linked_ID)
        #keep Count data only for each comparison
        needed<-which(colnames(df) %in% rownames(coldata.1))    
        df.1 <- df[,needed]
        needed2<-which(colnames(df) %in% rownames(coldata.2))    
        df.2 <- df[,needed2]
        #Convert Resuts table into a data.frame
        res2 <- as.data.frame(res2)
        #decide what otu to save 
        otu.to.save <-as.character(rownames(res2))
        #from relative table we should get the mean across the row of the otu table
        df.1.meanRA <- rowMeans(df.1)
        df.2.meanRA <- rowMeans(df.2)
        #need to subset AND reorder just the otus that we have 
        df.1.meanRA.save <- df.1.meanRA[otu.to.save]
        df.2.meanRA.save <- df.2.meanRA[otu.to.save]
        #add the abundnace data for the res dataframe
        res2$abundance.1 <- df.1.meanRA.save
        res2$abundance.2 <- df.2.meanRA.save
        #Set Names of Results Table
        res2 <- setNames(cbind(rownames(res2), res2, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
        #Get Assay Data For Compairson 3
        #GenusData <-as.data.frame(assay(ddsvbal3)) #pruned to selected Genuses based on abundance
        ##Create Relative Abundance Table
        #df <-
        #    GenusData %>% 
        #        rownames_to_column('gs') %>%
        #        group_by(gs) %>% 
        #        summarise_all(funs(sum)) %>%
        #        mutate_if(is.numeric, funs(./sum(.))) %>%
        #        column_to_rownames('gs')
        #Get the ColData for Each Comparison
        coldata.1 <- coldata2[coldata2$three_groups=="Cluster_2",] %>%
                    select(Study_Linked_ID)
        coldata.2 <- coldata2[coldata2$three_groups=="Cluster_3",] %>%
                    select(Study_Linked_ID)
        #keep Count data only for each comparison
        needed<-which(colnames(df) %in% rownames(coldata.1))    
        df.1 <- df[,needed]
        needed2<-which(colnames(df) %in% rownames(coldata.2))    
        df.2 <- df[,needed2]
        #Convert Resuts table into a data.frame
        res3 <- as.data.frame(res3)
        #decide what otu to save 
        otu.to.save <-as.character(rownames(res3))
        #from relative table we should get the mean across the row of the otu table
        df.1.meanRA <- rowMeans(df.1)
        df.2.meanRA <- rowMeans(df.2)
        #need to subset AND reorder just the otus that we have 
        df.1.meanRA.save <- df.1.meanRA[otu.to.save]
        df.2.meanRA.save <- df.2.meanRA[otu.to.save]
        #add the abundnace data for the res dataframe
        res3$abundance.1 <- df.1.meanRA.save
        res3$abundance.2 <- df.2.meanRA.save
        #Set Names of Results Table
        res3 <- setNames(cbind(rownames(res3), res3, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
        #Write Tables of Differential Analysis
        write.table(res1,file=paste0(counts,".Cluster_1_vs_Cluster_3_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        write.table(res2,file=paste0(counts,".Cluster_1_vs_Cluster_2_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        write.table(res3,file=paste0(counts,".Cluster_2_vs_Cluster_3_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        #Remove Any Data without LOGFC data                              
        res1 <- res1[!is.na(res1$logFC),]
        res2 <- res2[!is.na(res2$logFC),]
        res3 <- res3[!is.na(res3$logFC),]
        #Match the three tables with the same Taxa
        needed<-which(res1$Gene.symbol %in% res3$Gene.symbol)    
        res1 <- res1[needed,]
        needed<-which(res2$Gene.symbol %in% res1$Gene.symbol)  
        res2 <- res2[needed,]
        needed<-which(res3$Gene.symbol %in% res1$Gene.symbol)  
        res3 <- res3[needed,]
        needed<-which(res1$Gene.symbol %in% res2$Gene.symbol)    
        res1 <- res1[needed,]
        needed<-which(res3$Gene.symbol %in% res2$Gene.symbol)    
        res3 <- res3[needed,]
        # Reorder Results based on FDR for comparison 1
        res1 = res1[order(res1$adj.P.Val, na.last = TRUE), ]
        # Keep only top 10
        res1order <- res1 %>% slice(1:10)
        # Reorder Results based on FDR for comparison 2
        res2 = res2[order(res2$adj.P.Val, na.last = TRUE), ]
        # Keep only top 10
        res2order <- res2 %>% slice(1:10)
        # Reorder Results based on FDR for comparison 3
        res3 = res3[order(res3$adj.P.Val, na.last = TRUE), ]
        # Keep only top 10
        res3order <- res3 %>% slice(1:10)
        #Combine order of three tables
        resorder <- rbind(res1order,res2order,res3order)
        #Remove any duplicates
        resorder <- resorder[!duplicated(resorder$Gene.symbol), ]
        #Keep only matching Taxa for all three tables
        res1 <- res1[res1$Gene.symbol %in% resorder$Gene.symbol,]
        #Set The order
        res1 <- res1[ order(match(res1$Gene.symbol, resorder$Gene.symbol)), ]
        #Create Variable for order
        res1 <- res1 %>% mutate(start = 1:n())
        #Keep only matching Taxa for all three tables
        res2 <- res2[res2$Gene.symbol %in% resorder$Gene.symbol,]
        #Set The order
        res2 <- res2[ order(match(res2$Gene.symbol, resorder$Gene.symbol)), ]
        #Keep only matching Taxa for all three tables
        res2 <- res2 %>% mutate(start = 1:n())
        #Keep only matching Taxa for all three tables
        res3 <- res3[res3$Gene.symbol %in% resorder$Gene.symbol,]
        #Set The order
        res3 <- res3[ order(match(res3$Gene.symbol, resorder$Gene.symbol)), ]
        #Keep only matching Taxa for all three tables
        res3 <- res3 %>% mutate(start = 1:n())
        #Set Variable for the three comparisson
        res1$group <- "A"
        res2$group <- "B"
        res3$group <- "C"
        #Convert Important columns to Numeric
        res1$adj.P.Val <-   as.numeric(as.character(res1$adj.P.Val))
        res1$logFC <-       as.numeric(as.character(res1$logFC))
        res1$abundance.1 <- as.numeric(as.character(res1$abundance.1))
        res1$abundance.2 <- as.numeric(as.character(res1$abundance.2))
        res2$adj.P.Val <-   as.numeric(as.character(res2$adj.P.Val))
        res2$logFC <-       as.numeric(as.character(res2$logFC))
        res2$abundance.1 <- as.numeric(as.character(res2$abundance.1))
        res2$abundance.2 <- as.numeric(as.character(res2$abundance.2))
        res3$adj.P.Val <-   as.numeric(as.character(res3$adj.P.Val))
        res3$logFC <-       as.numeric(as.character(res3$logFC))
        res3$abundance.1 <- as.numeric(as.character(res3$abundance.1))
        res3$abundance.2 <- as.numeric(as.character(res3$abundance.2))
        #Bind the three Comparison TAbles
        resy <- rbind(res1,res2,res3)
        #Replace NA
        resy <- resy %>% mutate(adj.P.Val = if_else(is.na(adj.P.Val), 0.9, adj.P.Val))
        #Create Variable for Color based on Comparison, FDR and LOGFC
        resy$col <- ifelse(resy$group=="A" & resy$adj.P.Val<0.2 & resy$logFC>0, "A",
                    ifelse(resy$group=="A" & resy$adj.P.Val<0.2 & resy$logFC<0, "B",
                    ifelse(resy$group=="B" & resy$adj.P.Val<0.2 & resy$logFC>0, "C",
                    ifelse(resy$group=="B" & resy$adj.P.Val<0.2 & resy$logFC<0, "B",
                    ifelse(resy$group=="C" & resy$adj.P.Val<0.2 & resy$logFC>0, "A",
                    ifelse(resy$group=="C" & resy$adj.P.Val<0.2 & resy$logFC<0, "C","D"))))))
        #PLOT IT
            ggsave(filename=paste0(counts,".BAL_three_groups_DESEQ2.pdf"),
            ggplot(resy, aes(y=reorder(Gene.symbol,-start), x=logFC,fill=col)) +
            facet_grid(~ group, scales = "free_y")+
            #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 500 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 500 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Bacterial Size
            #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 5000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 5000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #DNA Bacterial Size
            geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 1000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 1000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Virome Size
            #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 1000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 1000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Fungi Size
            #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 10000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 10000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+
            #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 200000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 200000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Phage size
            #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 200000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 200000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #DNA Fungi Size
            #scale_fill_manual(values=c("#D01C8B","white"))+ #RNA Virome Colors
            scale_fill_manual(values=c("#D01C8B","#4DAC26","#FFD479","white"))+ #RNA Fungi Colors/Bacterial Colors
            #scale_fill_manual(values=c("#D01C8B","#4DAC26","white"))+ #RNA Fungi Colors/Bacterial Colors
            #scale_fill_manual(values=c("#D01C8B","#FFD479","white"))+ #RNA Phages Colors
            #scale_fill_manual(values=c("#FFD479","white"))+ 
            #scale_fill_manual(values=c("#D01C8B","white"))+ #DNA Fungi Colors
            #scale_fill_manual(values=c("white"))+ #DNA Virome Colors/DNA Phage Colors
            scale_size_continuous(range=c(1, 27),guide=FALSE)+
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
            width=20, height=5)
}