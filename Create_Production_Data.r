## THIS FUNCTION CREATES MRNA AND MIRNA DATASETS THAT ARE PROCESSED TO:
## 1. EXCLUDE NORMAL TISSUE SAMPLES
## 2. AVERAGE SAMPLES THAT ARE FROM THE SAME PHYSICAL SAMPLE
## 3. REMOVE COUNT INFORMATION FROM GENES/MIRNA THAT HAVE TOO FEW READS TO BE
## USEFUL TO ANALYZE
##
## THE RESULTING DATASETS ARE THEN PRINTED IN THE DIRECTORIES THAT THE ORIGINAL
## DATA ORIGINATED AND GIVEN AN _USE.TXT SUFFIX
library(edgeR)
library(peer)
library(mvoutlier)
library(ggplot2)
library(GGally)
library("biomaRt")
library(Homo.sapiens)
    
source("/group/stranger-lab/askol/TCGA/Code/Create_Production_Data_funcs.r")

Create_mRNA_Data <- function(project, use.existing.gene.info = TRUE){

    ## MUCH OF WHAT IS BELOW IS DERIVED FROM:
    ## https://www.bioconductor.org/help/workflows/RNAseq123/
    setwd("/group/stranger-lab/askol/TCGA/TCGA_Expression/")
    print(date())
    
  
    EXPDIR = "/group/stranger-lab/askol/TCGA/TCGA_Expression/"
    METADIR = "/group/stranger-lab/askol/TCGA/TCGA_MetaData/"
    
    ## GET GENE INFORMATION FROM BIOMART ##
    if (use.existing.gene.info==TRUE){
        gene.info <- read.table(file = "Gene.Info.txt", as.is=T, header=T, sep="\t")
        print("Using existing gene.info file")
    }else{
        gene.info <- get.gene.info()
        write.table(file = "Gene.Info.txt", gene.info, quote=F, row.names=F,
                    col.names=T, sep="\t")
    }
    
    ## GET INFORMATION ABOUT SAMPLES (CASE.ID, GENDER)
    map <- make.m.map(project)
    
    countFile <- paste0(EXPDIR,project, ".DATA.txt")
    fpkmFile <- paste0(EXPDIR,"FPKM-UQ/",project,"/",project,"-FPKM-UQ.txt")

    CountData <- get.m.data(countFile, fpkmFile, map, gene.info)

    ## REMOVE SAMPLES WITH MISSING VALUES OF THE COV ##
    ind.na <- which(is.na(map$cases.0.demographic.gender))
    ids.na <- map$cases.0.case_id[ind.na]

    print(paste0("Removing ",length(ids.na), " cases with no gender information"))
    
    if (length(ind.na) > 0){
        map <- map[-ind.na,]
        ind.na = which(gsub("A","", names(CountData)) %in% ids.na)
        CountData = CountData[,-ind.na]
    }
    map$case.id <- paste0("A",map$cases.0.case_id)
    gender <- map[,c("case.id", "gender")]
    rm.dupes <- which(duplicated(gender$case.id))
    if (length(rm.dupes) > 0){
        gender <- gender[-rm.dupes,]
    }
    group <- data.frame(case.id = names(CountData)[-c(1:2)])
    group <- merge(group, gender, by.x=1,
                   by.y="case.id", all.x = TRUE, all.y=FALSE, sort=F)
    rownames(CountData) <- CountData$ensmbl

    rm.genes <- which(rowSums(is.na(CountData[,-c(1:2)])) == (ncol(CountData)-2))
    
    ## GET GENE INFORMATION
    geneid <- rownames(CountData)
    genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                    keytype="ENSEMBL")
    xy.genes <- genes$ENSEMBL[genes$TXCHROM %in% c("chrX", "chrY")]
    xy.ind <- which(rownames(CountData) %in% xy.genes)
    rm.genes.xy <- unique(c(xy.ind, rm.genes))

    ## d.auto IS TO BE USED FOR AUTOSOMAL GENES (NORMALIZATION WAS DONE WITHOUT
    ## SEX CHROMOSOMES ##
    ## create the DGEList object and calculate variance
    d.auto <- DGEList(counts = CountData[-rm.genes.xy , -c(1:2)], group=group$gender)
    d.auto <- calcNormFactors(d.auto, method="TMM")
    d.auto <- estimateCommonDisp(d.auto,verbose=TRUE)
    d.auto <- estimateTagwiseDisp(d.auto)

    ## D.ALL IS TO BE USED FOR THE SEX CHROMOSOMES (NORMALIZATION WAS DONE WITH
    ## ALL GENES ON ALL CHROMOSOMES ##
    d.all <- DGEList(counts = CountData[-rm.genes , -c(1:2)], group=group$gender)
    d.all <- calcNormFactors(d.all, method="TMM")
    d.all <- estimateCommonDisp(d.all,verbose=TRUE)
    d.all<- estimateTagwiseDisp(d.all)  

    xy.ind <- rownames(d.all$counts) %in% xy.genes
    d.xy <- d.all[xy.ind, ,keep.lib.sizes=FALSE]

    ## CONVERT COUNTS TO COUNTS PER MISSION AND LOG COUNT PER MILLION
    cpm.auto <- cpm(d.auto)
    lcpm.auto <- cpm(d.auto, log=TRUE) ## NOTE: LOG=TRUE ADDS 0.25 TO EACH COUNT TO AVOID 0
    cpm.xy <- cpm(d.xy)
    lcpm.xy <- cpm(d.xy, log=TRUE) ## NOTE: LOG=TRUE ADDS 0.25 TO EACH COUNT TO AVOID 0

    
    ## DISTRIBUTION OF COUNT DATA ##
    plot.file <- paste0(project, "_mRNA_count_dist.pdf")
    pdf(file = plot.file)
    library(RColorBrewer)
    nsamples <- ncol(d.auto)
    col <- brewer.pal(nsamples, "Paired")
    plot(density(lcpm.auto[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
         main="", xlab="")
    title(main="Filtered data", xlab="Log-cpm")
    abline(v=0, lty=3)
    for (i in 2:nsamples){
        den <- density(lcpm.auto[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
    }
    dev.off()

    ## CALCULATE PEERS ##
    ## EXPECTING ROWS = SAMPLES, COLUMNS = GENES ##
    model <- PEER()
    PEER_setPhenoMean(model,as.matrix(t(lcpm.auto)))
    PEER_setCovariates(model, as.matrix(1*(group$gender=="male")))
                       
    ## SET 10 CONFOUNDERS ##
    PEER_setNk(model,10)

    ## RUN INFERENCE ##
    PEER_update(model)
                       
    factors = PEER_getX(model)
    weights = PEER_getW(model)
    precision = PEER_getAlpha(model)
    residuals = PEER_getResiduals(model)

    factors <- data.frame(factors[,-1]) ## remove sex since alread in d$samples
    names(factors) <- paste0("PEER",1:ncol(factors))
    ## ADD FIRST 10 PEERS TO D$SAMPLE
    d.auto$samples <- cbind(d.auto$samples, factors)
    d.xy$samples <- cbind(d.xy$samples, factors)
    
    ## TEST FOR ASSOCIATION BETWEEN PEER AND SEX ##
    dat <- data.frame(factors, sex=d.auto$sample$group)
    names(dat)[ncol(dat)] = "sex"
    lms <- lapply(grep("PEER", names(dat)),
                  function(x) as.vector(t(
                          summary(lm(dat[,x] ~ dat$sex))$coef[2,])))
    lms.mtx <- do.call(rbind, lms)

    ## REMOVE PEER IF P<0.01 ##
    rm.peer <- which(lms.mtx[,4] < 0.01)
    if (length(rm.peer) >0){
        d.auto$samples <- d.auto$samples[,-which(names(d.auto$samples)
                                                 %in% paste0("PEER",rm.peer))]
        d.xy$samples <- d.xy$samples[,-which(names(d.xy$samples)
                                               %in% paste0("PEER",rm.peer))]
    }
    
    ## SAVE DATA ##
    file <- paste0("Limma_Voom_Use/",project,"_mrna_autosome_DGEList.RDS")
    saveRDS(d.auto, file=file)
    file <- paste0("Limma_Voom_Use/",project,"_mrna_sex_DGEList.RDS")
    saveRDS(d.xy, file=file)

    print(paste0("Wrote ", file))
    print(paste0("Wrote ", gsub("sex","autosome", file)))
    
}

Create_miRNA_Data <- function(project){
        
    ROOT <- "/group/stranger-lab/askol/TCGA/"
    MIRDIR <- "/group/stranger-lab/askol/TCGA/TCGA_microRNA/"
    setwd(MIRDIR)

    map <- make.mi.map(project)
    
    mi.file <- paste0(MIRDIR,project,"/",project,"_miR_count.txt")
    
    mi <- get.mi.data(mi.file, map)
    mi$ID <- gsub("-","\\.",mi$ID)
    rownames(mi) <- mi$ID
    
    ## REMOVE SAMPLES WITH MISSING VALUES OF THE COV ##
    ind.na <- which(is.na(map$cases.0.demographic.gender))
    ids.na <- map$cases.0.case_id[ind.na]

    print(paste0("Removing ",length(ids.na), " cases with no gender information"))
    
    if (length(ind.na) > 0){
        map <- map[-ind.na,]
        ind.na = which(gsub("A","", names(mi)) %in% ids.na)
        mi = mi[,-ind.na]
    }

    ## ### --------- ######
    map$case.id <- paste0("A",map$cases.0.case_id)
    
    group <- data.frame(case.id = names(mi)[-1])
    group <- merge(group, map[,c("case.id", "gender")],by.x=1,
                   by.y="case.id", sort=F)
    rownames(mi) <- mi$ID

    ## REMOVE MIRNAS WHERE ALL VALUES ARE NA
    rm.genes <- which(rowSums(is.na(mi[,-1])) == (ncol(mi)-1))
    mi <- mi[-rm.genes,]
    xy.mirs <- get.xy.mirs()

    ## SPLIT MI INTO AUTOSOME AND XY ##
    mi.xy <- mi[mi$ID %in% xy.mirs, -1]
    mi.auto <- mi[!mi$ID %in% xy.mirs, -1]

    ## d.auto IS TO BE USED FOR AUTOSOMAL GENES (NORMALIZATION WAS DONE WITHOUT
    ## SEX CHROMOSOMES ##
    ## create the DGEList object and calculate variance
    d.auto <- DGEList(counts = mi.auto, group=group$gender)
    d.auto <- calcNormFactors(d.auto, method="TMM")
    d.auto <- estimateCommonDisp(d.auto,verbose=TRUE)
    d.auto <- estimateTagwiseDisp(d.auto)

    ## D.ALL IS TO BE USED FOR THE SEX CHROMOSOMES (NORMALIZATION WAS DONE WITH
    ## ALL GENES ON ALL CHROMOSOMES ##
    d.all <- DGEList(counts = mi[-rm.genes , -1], group=group$gender)
    d.all <- calcNormFactors(d.all, method="TMM")
    d.all <- estimateCommonDisp(d.all,verbose=TRUE)
    d.all<- estimateTagwiseDisp(d.all)  

    xy.ind <- rownames(d.all$counts) %in% xy.mirs
    d.xy <- d.all[xy.ind, ,keep.lib.sizes=FALSE]

    ## CONVERT COUNTS TO COUNTS PER MISSION AND LOG COUNT PER MILLION
    cpm.auto <- cpm(d.auto)
    lcpm.auto <- cpm(d.auto, log=TRUE) ## NOTE: LOG=TRUE ADDS 0.25 TO EACH COUNT TO AVOID 0
    cpm.xy <- cpm(d.xy)
    lcpm.xy <- cpm(d.xy, log=TRUE) ## NOTE: LOG=TRUE ADDS 0.25 TO EACH COUNT TO AVOID 0

    
    ## DISTRIBUTION OF COUNT DATA ##
    plot.file <- paste0(project, "_mir_count_dist.pdf")
    pdf(file = plot.file)
    library(RColorBrewer)
    nsamples <- ncol(d.auto)
    col <- brewer.pal(nsamples, "Paired")
    plot(density(lcpm.auto[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
         main="", xlab="")
    title(main="Filtered data", xlab="Log-cpm")
    abline(v=0, lty=3)
    for (i in 2:nsamples){
        den <- density(lcpm.auto[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
    }
    dev.off()

    ## CALCULATE PEERS ##
    ## EXPECTING ROWS = SAMPLES, COLUMNS = GENES ##
    model <- PEER()
    PEER_setPhenoMean(model,as.matrix(t(lcpm.auto)))
    PEER_setCovariates(model, as.matrix(1*(group$gender=="male")))
                       
    ## SET 10 CONFOUNDERS ##
    PEER_setNk(model,10)

    ## RUN INFERENCE ##
    PEER_update(model)
                       
    factors = PEER_getX(model)
    weights = PEER_getW(model)
    precision = PEER_getAlpha(model)
    residuals = PEER_getResiduals(model)

    factors <- data.frame(factors[,-1]) ## remove sex since alread in d$samples
    names(factors) <- paste0("PEER",1:ncol(factors))
    ## ADD FIRST 10 PEERS TO D$SAMPLE
    d.auto$samples <- cbind(d.auto$samples, factors)
    d.xy$samples <- cbind(d.xy$samples, factors)
    
    ## TEST FOR ASSOCIATION BETWEEN PEER AND SEX ##
    dat <- data.frame(factors, sex=d.auto$sample$group)
    names(dat)[ncol(dat)] = "sex"
    lms <- lapply(grep("PEER", names(dat)),
                  function(x) as.vector(t(
                          summary(lm(dat[,x] ~ dat$sex))$coef[2,])))
    lms.mtx <- do.call(rbind, lms)

    ## REMOVE PEER IF P<0.01 ##
    rm.peer <- which(lms.mtx[,4] < 0.01)
    if (length(rm.peer) >0){
        d.auto$samples <- d.auto$samples[,-which(names(d.auto$samples)
                                                 %in% paste0("PEER",rm.peer))]
        d.xy$samples <- d.xy$samples[,-which(names(d.xy$samples)
                                               %in% paste0("PEER",rm.peer))]
    }
    
    ## SAVE DATA ##
    file <- paste0("Limma_Voom_Use/",project,"_mirna_autosome_DGEList.RDS")
    saveRDS(d.auto, file=file)
    file <- paste0("Limma_Voom_Use/",project,"_mirna_sex_DGEList.RDS")
    saveRDS(d.xy, file=file)

    print(paste0("Wrote ", file))
    print(paste0("Wrote ", gsub("sex","autosome", file)))

}

Create_Production_Data <- function(project){
    
    print(paste0("Creating mRNA file for project ",project))
    Create_mRNA_Data(project)
    
    print(paste0("Creating miRNA file for project ",project))
    Create_miRNA_Data(project)
    
    print(paste0("Done with project ",project))

}

## ------ MAIN --------##
args <- commandArgs(TRUE)
project = args[1]
Create_Production_Data(project)
print(date())

finish.file = paste0(project,".finished")
system(paste0("touch ",finish.file))
