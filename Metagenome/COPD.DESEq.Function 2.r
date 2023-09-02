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
        coldata$Study_Linked_ID <- coldata$GTC_UPDATED_NAME_MetagenomeMH2021 #Metagenome
        #coldata$Study_Linked_ID <- coldata$GTC_Updated_Name_MetatranscriptomeMH2021 #Metatranscriptome
        #coldata$Study_Linked_ID <- gsub("-","_",coldata$Study_Linked_ID) #Metatranscriptome
        #coldata$Study_Linked_ID <- gsub(".Sup","_Sup",coldata$Study_Linked_ID) #Metatranscriptome
        coldata$Study_Linked_ID <- gsub("-",".",coldata$Study_Linked_ID) #Metagenome
        coldata <- coldata[!duplicated(coldata$Study_Linked_ID), ]
        coldata <- coldata[coldata$Study_Linked_ID!="n.a", ]
        rownames(coldata) <- coldata$Study_Linked_ID
        coldata$Sample.Type <- coldata$Sample_type_IS
        coldata$Sample.Type <- ifelse(coldata$Sample.Type=="BKG","UA",
                        ifelse(coldata$Sample.Type=="UA","BKG",as.character(coldata$Sample.Type))) #Metagenome
        #Load Count Data
        mycounts  <-read.delim2(paste0(counts,".txt"), sep="\t", row.names=1)
        #mycounts   <- mycounts %>% select(!category) #BacterialData      
        #Create Relative Abundance Table
        rel <-as.data.frame(mycounts) #make counts dataframe
        rel = data.frame(lapply(rel, function(x) as.numeric(as.character(x))),
                           check.names=F, row.names = rownames(rel))
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
        #Remove 0 Counts
        i <- (colSums(d1, na.rm=T) != 0) # T if colSum is not 0, F otherwise
        d1 <- d1[ , i]
        wanted3<-which(rownames(coldata2) %in% colnames(d1))
        coldata2<-coldata2[wanted3,]
        #DESEQ Object
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
        #Transforming data - Option #2 is variance stabilizing transformation
        vsdv <- varianceStabilizingTransformation(ddsv)
        #vsdv <- rlog(ddsv, fitType="local",blind=FALSE)
        #----------------------
        ##HEATMAP
        #----------------------
        #filter out genes where there are less than 3 samples with normalized counts greater than or equal to 100.
        idx <- rowSums(relcounts3 >= 0.0003 ) >= 9
        df <- relcounts3[idx,]
        heatmap_raw <- mycounts3[idx,]
        #cluster Genuses(row)
        GenusData.Bray.dist <-vegdist(df, method = "bray")
        #cluster samples(Col)
        Samples.Bray.dist = vegdist(t(df), method="bray")
        #Set Color Scale for Heatmap
        mypalette <- colorRampPalette(c('#ffffff','#4169E1','#0000CD'))(100)
        #Set Colors for each sample type for HeatMap
        annon_colors= list(Sample.Type=c(BKG="#BEBEBE", BAL="#1700F5", UA="#F7A501"))
        #Choose lables for Samples
        df2 <- data.frame(Sample.Type = colData(ddsv)[,c("Sample.Type")], row.names = rownames(colData(ddsv)))
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
            height = 15, width = 20)
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
        Reads$or <-ifelse(Reads$Sample.Type=="BKG", 1,NA)
        Reads$or <-ifelse(Reads$Sample.Type=="BAL",2 ,Reads$or)
        Reads$or <-ifelse(Reads$Sample.Type=="UA",3 ,Reads$or)
        #Create Figure
            ggsave(filename=paste0(counts,".Reads.pdf"),
            ggplot(Reads, aes(x= reorder(Sample.Type, +or), y=Reads, fill=Sample.Type)) + 
            stat_boxplot(geom ='errorbar', width=0.1)+
            geom_boxplot(outlier.shape = NA, width=0.5)+
            geom_jitter(shape=1, position=position_jitter(0.2))+
            scale_fill_manual(values=c("#1700F5", "#BEBEBE","#F7A501")) +
            scale_x_discrete(labels = c('BKG','BAL','UA'))+ 
            scale_y_continuous(name="Reads",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
            xlab("Sample Type")+
            theme,
            height = 7, width = 5)

        #Keep only BAL
        bal <- coldata2[coldata2$Sample.Type=="BAL",]
        #keep only matching IDs from count data
        wanted<-which(colnames(heatmap_raw) %in% rownames(bal))
        wanted2<-which(colnames(df) %in% rownames(bal))
        #keep only matching IDs from count data
        heatmap_raw3<-heatmap_raw[,wanted]
        heatmap3<-df[,wanted2]
        write.table(heatmap_raw3,file="Metagenome_BAL_raw_pruned_count_table.txt", sep="\t", col.names = NA, row.names = TRUE)
        write.table(heatmap3,file="Metagenome_BAL_Rel_Abundance_pruned_count_table.txt", sep="\t", col.names = NA, row.names = TRUE)

        #BAL Heatmap
        #cluster Genuses(row)
        GenusData.Bray.dist <-vegdist(heatmap3, method = "bray")
        #cluster samples(Col)
        Samples.Bray.dist = vegdist(t(heatmap3), method="bray")
        #Set Color Scale for Heatmap
        mypalette <- colorRampPalette(c('#ffffff','#4169E1','#0000CD'))(100)
        #Set Colors for each sample type for HeatMap
        annon_colors= list(Subject_Type = c(Smoker.Control="#4DAC26", COPD="#D01C8B"))
        #Choose lables for Samples
        df2 <- data.frame(Subject_Type = coldata2$Subject_Type, row.names = coldata2$Study_Linked_ID)
        #convert to data.frame
        #Create dataframe of count data
        heatmap3 <- as.data.frame(heatmap3)
        #Create Heatmap
        pdf("Metagenome_pruned_BAL_Heatmap.pdf", height = 15, width = 20)
            pheatmap(heatmap3, cluster_rows=TRUE, show_rownames=TRUE, 
            cluster_cols=TRUE,annotation_col=df2,
            scale="row",
            clustering_distance_rows = GenusData.Bray.dist,
            clustering_distance_cols = Samples.Bray.dist,
            clustering_method="mcquitty",
            gaps_col=50,
            border_color="black",
            color = colorRampPalette(c('#4169E1','#ffffff','#0000CD'))(100),
            annotation_colors=annon_colors,legend=FALSE)
        dev.off()
        #----------------------
        ##PCOA PLOT
        #----------------------
        #Create Distance Matrix
        vsdv0 <- ifelse(assay(vsdv)<0,0,assay(vsdv))
        vegdist   = vegdist(t(vsdv0), method="bray")
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
        #Run the Statistics
        samplepermanova <- adonis(vegdist ~ Sample.Type, data.adonis)
        samplepermanova <- as.data.frame(samplepermanova$aov.tab)
        samplepermanova <- samplepermanova$'Pr(>F)'[1]
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
        ##Alpha Diversity
        #----------------------
        #Calcultes Shannon Diversity
        ddsv$Shannon = diversity(vsdv0, index = "shannon", MARGIN = 2, base = exp(1))
        #ddsv$Shannon = diversity(d2, index = "shannon", MARGIN = 2, base = exp(1))
        #Convert to data frame for ggplot
        shannon = as.data.frame(colData(ddsv))
        #Remove any zero values
        shannon[shannon==0] <- NA
        #Set Order Of Figure
        shannon$or <-ifelse(shannon$Sample.Type=="BKG", 1,NA)
        shannon$or <-ifelse(shannon$Sample.Type=="BAL",2 ,shannon$or)
        shannon$or <-ifelse(shannon$Sample.Type=="UA",3 ,shannon$or)
        #Check Statistics
        sampleshannon <- kruskal.test(Shannon ~ Sample.Type, data = shannon)
        sampleshannon <- sampleshannon$p.value
        #Plot It
            ggsave(filename=paste0(counts,".SHANNON_vsd_sampletype_Pvalue_",sampleshannon,".pdf"),
            ggplot(shannon, aes(x= reorder(Sample.Type, +or), y=Shannon, fill=Sample.Type)) + 
            stat_boxplot(geom ='errorbar', width=0.1)+
            geom_boxplot(outlier.shape = NA, width=0.5)+
            #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
            geom_jitter(shape=1, position=position_jitter(0.2))+
            scale_fill_manual(values=c("#1700F5", "#BEBEBE", "#F7A501")) +
            scale_x_discrete(labels = c('BKG','BAL','UA'))+ 
            #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
            #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
            #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
            ylab("Shannon Diversity") + 
            xlab("Sample Type")+
            #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
            #geom_point(color=cols) +
            #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
            theme, 
            height = 7, width = 5)
        #----------------------
        ##DECONTAM
        #----------------------
        GenusData <-as.data.frame(assay(ddsv)) #pruned to selected Genuses based on abundance        
        #Create Relative Abundance Table
        df <-
            GenusData %>% 
                rownames_to_column('gs') %>%
                group_by(gs) %>% 
                summarise_all(funs(sum)) %>%
                mutate_if(is.numeric, funs(./sum(.))) %>%
                column_to_rownames('gs')
        #Create a TRUE/FALSE Variable for Subject and Control
        coldata2$is.neg <- coldata2$Sample.Type=="BKG"
        #Transpose to matrix
        otus <- t(as.matrix(df))
        ## extract only those samples in common between the two tables
        common.sample.ids <- intersect(rownames(coldata2), rownames(otus))
        otus <- otus[common.sample.ids,]
        coldata3 <- coldata2[common.sample.ids,]
        #Check Contaminants based on BKG samples
        contamdf.prev <- isContaminant(otus, method="prevalence", neg=coldata3$is.neg,threshold=0.5)
        table(contamdf.prev$contaminant)
        #Select only BKG Samples
        BKG <- df %>% select(rownames(coldata3[coldata3$Sample.Type=='BKG',]))
        #List of the Contaminants Identified
        otu.to.save <-as.character(rownames(contamdf.prev[contamdf.prev$contaminant==TRUE,]))
        #sum of BKG Abundances
        BKG <-rowSums(BKG)
        #Only select contaminants
        BKG <- BKG[otu.to.save]
        #Sort table by Total Abundance
        BKG <- sort(BKG,decreasing=TRUE)
        #transform to Data Frame
        BKG <- as.data.frame(BKG)
        #Select top 10
        BKG <- BKG %>% dplyr::slice(1:10)
        #write Table
        write.table(BKG,file=paste0(counts,".top10_contaminants.txt"), sep="\t", col.names = NA, row.names = TRUE)
        #----------------------
        ##DESEQ
        #----------------------
        #Create Deseq object for BAL analysis
        ddsvbal <- DESeqDataSetFromMatrix(countData = d1,
                              colData = coldata2,
                              design= ~ Sample.Type)
        #Subset just the BAL
        ddsvbal <- ddsvbal[, ddsvbal$Sample.Type %in% "BAL"]
        #Subset Rel Table for BAL only
        wanted<-which(colnames(d2) %in% colnames(ddsvbal))
        #keep only matching IDs from count data
        d3<-d2[,wanted]
        #Covert Variable to Factor
        ddsvbal$Subject_Type <- as.factor(ddsvbal$Subject_Type)
        #Create Seperate tables with each of the different compairsons
        #ddsvbal1 <- ddsvbal[, ddsvbal$three_groups %in% c("Dead","Less_Than_28_days_on_vent")]
        #ddsvbal2 <- ddsvbal[, ddsvbal$three_groups %in% c("Less_Than_28_days_on_vent","Greater_Than_28_days_on_vent")]
        #ddsvbal3 <- ddsvbal[, ddsvbal$three_groups %in% c("Greater_Than_28_days_on_vent","Dead")]
        #Create Tables of First Comparison
        balcounts <- assay(ddsvbal)
        balmeta   <- colData(ddsvbal)
        data_keep_rows <- c("Streptococcus mitis", "Veillonella parvula", "Prevotella melaninogenica")
        balcounts_subset <- balcounts[rownames(balcounts) %in% data_keep_rows,]
        transform <- t(balcounts_subset)
        colnames(transform) <- c("Streptococcus_mitis", "Veillonella_parvula", "Prevotella_melaninogenica")
        merged <- merge(balmeta,transform, by=0, all=TRUE)
        merged <- merged %>% select(Study_Linked_ID,Subject_Type,Streptococcus_mitis,Veillonella_parvula,Prevotella_melaninogenica)
        merged2 <- melt(merged, id.vars=c("Study_Linked_ID", "Subject_Type"))
        merged2$Subject_Type <- factor(merged2$Subject_Type, levels=c("Smoker.Control", "COPD"))
        merged2$variable <- factor(merged2$variable, levels=c("Streptococcus_mitis", "Veillonella_parvula", "Prevotella_melaninogenica"))

        #PLOT IT
        pdf("Metagenome_BAL_MOC.pdf", height = 10, width = 15)
            ggplot(merged2, aes(x= Subject_Type, y=value, fill=Subject_Type)) + 
            stat_boxplot(geom ='errorbar', width=0.1)+
            geom_boxplot(outlier.shape = NA, width=0.5)+
            geom_jitter(shape=1, position=position_jitter(0.2))+
            scale_fill_manual(values=c("#4DAC26","#D01C8B")) +
            facet_grid(cols=vars(variable)) + 
            #scale_x_discrete(labels = c('<28 Days','>28 Days','Deceased'))+ 
            #ylab("Shannon Diversity") + 
            xlab("Subject Type")+
            scale_y_continuous(name="Reads",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
            #scale_y_continuous(breaks = seq(2, 5, .5), limits = c(2, 5)) + #For Fungi
            #scale_y_continuous(breaks = seq(3, 5, .5), limits = c(3, 5)) + #For Fungi
            theme
        dev.off()
        #Check Statistics
        merged2_s <- merged2[merged2$variable=="Streptococcus_mitis",]
        merged2_v <- merged2[merged2$variable=="Veillonella_parvula",]
        merged2_p <- merged2[merged2$variable=="Prevotella_melaninogenica",]
        kruskal.test(value ~ Subject_Type, data = merged2_s)
        kruskal.test(value ~ Subject_Type, data = merged2_v)
        kruskal.test(value ~ Subject_Type, data = merged2_p)
        

        #write.table(balcounts,file="metagenome_bacteria_BAL.txt", sep="\t", col.names = NA, row.names = TRUE)
        ##Create Tables of First Comparison
        #balcounts1 <- assay(ddsvbal1)
        #balmeta1   <- colData(ddsvbal1)
        ##Create Tables of Second Comparison
        #balcounts2 <- assay(ddsvbal2)
        #balmeta2   <- colData(ddsvbal2)
        ##Create Tables of Third Comparison
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

        #Loadings
        #Set Variables
        counts <- balcounts
        groups <- as.data.frame(balmeta)
        ### Set taxa cutoffs for loadings
        n_samples <- nrow(groups)
        min_reads <- 10
        ###
        taxa <- rownames(counts)[rowSums(counts >= min_reads) >= n_samples]
        wascores <- wascores(cmdScale, as.data.frame(t(vsdv0[taxa,])))
        wascores <- wascores(CmdScale, as.data.frame(t(vsdvbal0)))
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
    pdf("Metagenome.BAL_PCOA_Loadings.pdf", height = 5, width = 10)
        plot_grid(plot.wascores.comp1, plot.wascores.comp2, labels = c('A', 'B'))    
    dev.off()
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
        df <- d3
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
            geom_point(color=cols, size = ifelse(res$logFC>=1 & res$adj.P.Val < alpha, 4000 * res$abundance.2, ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, 4000 * res$abundance.1,1)),alpha=0.6) + #Chose Colors for dots
            geom_text_repel(aes(label=ifelse(res$logFC<(-1) & res$adj.P.Val < alpha , as.character(res$Gene.symbol),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
            geom_text_repel(aes(label=ifelse(res$logFC>3 & res$adj.P.Val < alpha , as.character(res$Gene.symbol),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
            theme(legend.position = "none") +
            geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
            xlab("Effect size: log2(fold-change)") +
            ylab("-log10(adjusted p-value)") + 
            #ylim(0,20)+
            theme,
            width=5, height=5)
        
        
        
        
        #Remove Any Data without LOGFC data
        #es1 <- res1[!is.na(res1$logFC),]
        #es2 <- res2[!is.na(res2$logFC),]
        #es3 <- res3[!is.na(res3$logFC),]
        #Match the three tables with the same Taxa
        #eeded<-which(res1$Gene.symbol %in% res3$Gene.symbol)    
        #es1 <- res1[needed,]
        #eeded<-which(res2$Gene.symbol %in% res1$Gene.symbol)  
        #es2 <- res2[needed,]
        #eeded<-which(res3$Gene.symbol %in% res1$Gene.symbol)  
        #es3 <- res3[needed,]
        #eeded<-which(res1$Gene.symbol %in% res2$Gene.symbol)    
        #es1 <- res1[needed,]
        #eeded<-which(res3$Gene.symbol %in% res2$Gene.symbol)    
        #es3 <- res3[needed,]
        # Reorder Results based on FDR for comparison 1
        #es1 = res1[order(res1$adj.P.Val, na.last = TRUE), ]
        # Keep only top 10
        #es1order <- res1 %>% slice(1:10)
        # Reorder Results based on FDR for comparison 2
        #es2 = res2[order(res2$adj.P.Val, na.last = TRUE), ]
        # Keep only top 10
        #es2order <- res2 %>% slice(1:10)
        # Reorder Results based on FDR for comparison 3
        #es3 = res3[order(res3$adj.P.Val, na.last = TRUE), ]
        # Keep only top 10
        #es3order <- res3 %>% slice(1:10)
        #Combine order of three tables
        #esorder <- rbind(res1order,res2order,res3order)
        #Remove any duplicates
        #esorder <- resorder[!duplicated(resorder$Gene.symbol), ]
        #Keep only matching Taxa for all three tables
        #es1 <- res1[res1$Gene.symbol %in% resorder$Gene.symbol,]
        #Set The order
        #es1 <- res1[ order(match(res1$Gene.symbol, resorder$Gene.symbol)), ]
        #Create Variable for order
        #es1 <- res1 %>% mutate(start = 1:n())
        #Keep only matching Taxa for all three tables
        #es2 <- res2[res2$Gene.symbol %in% resorder$Gene.symbol,]
        #Set The order
        #es2 <- res2[ order(match(res2$Gene.symbol, resorder$Gene.symbol)), ]
        #Keep only matching Taxa for all three tables
        #es2 <- res2 %>% mutate(start = 1:n())
        #Keep only matching Taxa for all three tables
        #es3 <- res3[res3$Gene.symbol %in% resorder$Gene.symbol,]
        #Set The order
        #es3 <- res3[ order(match(res3$Gene.symbol, resorder$Gene.symbol)), ]
        #Keep only matching Taxa for all three tables
        #es3 <- res3 %>% mutate(start = 1:n())
        #Set Variable for the three comparisson
        #es1$group <- "A"
        #es2$group <- "B"
        #es3$group <- "C"
        #Convert Important columns to Numeric
        #es1$adj.P.Val <-   as.numeric(as.character(res1$adj.P.Val))
        #es1$logFC <-       as.numeric(as.character(res1$logFC))
        #es1$abundance.1 <- as.numeric(as.character(res1$abundance.1))
        #es1$abundance.2 <- as.numeric(as.character(res1$abundance.2))
        #es2$adj.P.Val <-   as.numeric(as.character(res2$adj.P.Val))
        #es2$logFC <-       as.numeric(as.character(res2$logFC))
        #es2$abundance.1 <- as.numeric(as.character(res2$abundance.1))
        #es2$abundance.2 <- as.numeric(as.character(res2$abundance.2))
        #es3$adj.P.Val <-   as.numeric(as.character(res3$adj.P.Val))
        #es3$logFC <-       as.numeric(as.character(res3$logFC))
        #es3$abundance.1 <- as.numeric(as.character(res3$abundance.1))
        #es3$abundance.2 <- as.numeric(as.character(res3$abundance.2))
        #Bind the three Comparison TAbles
        #esy <- rbind(res1,res2,res3)
        #Replace NA
        #esy <- resy %>% mutate(adj.P.Val = if_else(is.na(adj.P.Val), 0.9, adj.P.Val))
        #Create Variable for Color based on Comparison, FDR and LOGFC
        #esy$col <- ifelse(resy$group=="A" & resy$adj.P.Val<0.2 & resy$logFC>0, "A",
        #           ifelse(resy$group=="A" & resy$adj.P.Val<0.2 & resy$logFC<0, "B",
        #           ifelse(resy$group=="B" & resy$adj.P.Val<0.2 & resy$logFC>0, "C",
        #           ifelse(resy$group=="B" & resy$adj.P.Val<0.2 & resy$logFC<0, "B",
        #           ifelse(resy$group=="C" & resy$adj.P.Val<0.2 & resy$logFC>0, "A",
        #           ifelse(resy$group=="C" & resy$adj.P.Val<0.2 & resy$logFC<0, "C","D"))))))
        #PLOT IT
        #   ggsave(filename=paste0(counts,".BAL_three_groups_DESEQ2.pdf"),
        #   ggplot(resy, aes(y=reorder(Gene.symbol,-start), x=logFC,fill=col)) +
        #   facet_grid(~ group, scales = "free_y")+
        #   #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 500 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 500 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Bacterial Size
        #   #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 5000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 5000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #DNA Bacterial Size
        #   #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 1000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 1000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Virome Size
        #   #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 1000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 1000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Fungi Size
        #   #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 10000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 10000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+
        #   geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 200000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 200000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #RNA Phage size
        #   #geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 200000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 200000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+ #DNA Fungi Size
        #   #scale_fill_manual(values=c("#D01C8B","white"))+ #RNA Virome Colors
        #   #scale_fill_manual(values=c("#D01C8B","#4DAC26","#FFD479","white"))+ #RNA Fungi Colors/Bacterial Colors
        #   scale_fill_manual(values=c("#D01C8B","#FFD479","white"))+ #RNA Phages Colors
        #   #scale_fill_manual(values=c("#FFD479","white"))+ 
        #   #scale_fill_manual(values=c("#D01C8B","white"))+ #DNA Fungi Colors
        #   #scale_fill_manual(values=c("white"))+ #DNA Virome Colors/DNA Phage Colors
        #   scale_size_continuous(range=c(1, 27),guide=FALSE)+
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

}
