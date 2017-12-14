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
## library(Homo.sapiens)
    
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
        ##         gene.info <- get.gene.info()
        gene.info <- Get_Gene_Info()
        write.table(file = "Gene.Info.txt", gene.info, quote=F, row.names=F,
                    col.names=T, sep="\t")
    }

    ## GET GENE ANNOTATION INFORMATION ##
    gene.info <- Get_Gene_Info()
    
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
    group <- data.frame(case.id = names(CountData)[-1])
    group <- merge(group, gender, by.x=1,
                   by.y="case.id", all.x = TRUE, all.y=FALSE, sort=F)
    rownames(CountData) <- CountData$ensmbl 

    ## ONLY USE COUNTDATA FOR GENES WITH ANNOTATION INFORMATION ##
    ## MOST GENES NOT IN ANNOTATION THAT ARE IN COUNTDATA ARE LINCRNA AND OTHER
    ## SMALL AND LESS INTERESTING/NORMAL BEHAVING GENES
    ind <- which(CountData$ensmbl %in% gene.info$ENSEMBL)
    print(paste0( sum(unique(CountData$ensmbl) %in% gene.info$ENSEMBL == F),
                  " ENSEMBL IDs in count data, but not annotation. Genes removed. . ."))
    CountData <- CountData[ind,]
    ind <- which(gene.info$ENSEMBL %in% CountData$ensmbl)
    gene.info <- gene.info[ind,]

    ## DETERMINE GENES ON X/Y CHROMSOMES
    xy.genes <- gene.info$ENSEMBL[grep("X|Y", gene.info$chr)]
    xy.count.ind <- which(rownames(CountData) %in% xy.genes)
    xy.info.ind <- which(gene.info$ENSEMBL %in% xy.genes)
     
    ## d.auto IS TO BE USED FOR AUTOSOMAL GENES (NORMALIZATION WAS DONE WITHOUT
    ## SEX CHROMOSOMES ##
    ## create the DGEList object and calculate variance
    d.auto <- DGEList(counts = CountData[-xy.count.ind , -1], group=group$gender,
                      genes = gene.info[-xy.info.ind,])
    d.auto <- calcNormFactors(d.auto, method="TMM")
    d.auto <- estimateCommonDisp(d.auto,verbose=TRUE)
    d.auto <- estimateTagwiseDisp(d.auto)

    ## D.ALL IS TO BE USED FOR THE SEX CHROMOSOMES (NORMALIZATION WAS DONE WITH
    ## ALL GENES ON ALL CHROMOSOMES ##
    d.all <- DGEList(counts = CountData[, -1], group=group$gender)
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
    d.all <- DGEList(counts = mi[ , -1], group=group$gender)
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

Get_Gene_Info <- function(){
  
    genes <- get.gene.annotation()

    ## ADD INFORMATION ABOUT GENES THAT HAVE WEAK/STRONG INTERACTION WITH MIRNA ##
    mi.gene <- read.table("/group/stranger-lab/askol/TCGA/hsa_MTI-4.txt",
                          as.is=T, header=T, sep="\t")
    ## REMOVE 
    mi.gene$miRNA <- gsub("-3p|-5p","",mi.gene$miRNA)
    mi.gene$miRNA <- gsub("miR", "mir", mi.gene$miRNA)

   
    ## ADD INFORMATION ABOUT LOCATION TO MI.GENE ##
    mir.info <- get.mir.info()
    mi.gene <- merge(mi.gene, mir.info, by.x = "miRNA", by.y="mir",
                     all.x=T, all.y=F)

    mir.info.mirs <- strsplit(split="-", mir.info$mir)
    mir.info.mirs <- sapply(mir.info.mirs, function(x){paste(x[1:3], collapse="-")})
    ind <- which((mir.info$mir %in% mi.gene$miRNA) == F)
    sum(mir.info$mir %in% mi.gene$miRNA)
    
    
    ## CREATE STRONG AND WEAK MIR INTERACTION INDICATORS ON A GENE/MIR BASIS
    weak.ind <- grep("Weak", mi.gene$Support.Type)
    mi.gene$mir.weak <- mi.gene$mir.strong <- 0
    mi.gene$mir.weak[weak.ind] <- 1
    mi.gene$mir.strong[-weak.ind] <- 1

    ## CREATE A LIST OF WEAK AND STRONG INTERACTION GENES
    ## NOTE: WEAK GENES IS ANY GENE THAT DOESN'T HAVE A STRONG INDICATOR
    ## 
    strong.genes <- unique(mi.gene$Target.Gene[-weak.ind])
    weak.genes <- setdiff(unique(mi.gene$Target.Gene[weak.ind]),
                          strong.genes)
    genes$intx <- NA
    genes$intx[genes$SYMBOL %in% strong.genes] <- "strong"
    genes$intx[genes$SYMBOL %in% weak.genes] <- "weak"
    genes$intx[genes$SYMBOL %in% c(strong.genes, weak.genes) == FALSE] <- "none"  

    ## CREATE STRONG AND WEAK GENE-MIR INTACTION INDICATORS ##
    ## STRONG INTERACTIN GENES HAVE AT LEAST ONE STRONG INTERACTION
    mi.gene$gene.weak <- mi.gene$gene.strong <- 0
    mi.gene$gene.weak[mi.gene$Target.Gene %in% weak.genes] <- 1
    mi.gene$gene.strong[mi.gene$Target.Gene %in% strong.genes] <- 1

    ## ADD NAMES OF STRONG AND WEAK INTERACTING MIRS FOR EACH PROTEIN MIRS,
    ## THEIR CHROMOSOMES AND START POSITIONS
    ind.weak <- which(mi.gene$mir.weak == 1)
    ind.strong <- which(mi.gene$mir.strong == 1)
    mirs.weak <- with(mi.gene[ind.weak,],
                      tapply(miRNA, Target.Gene, paste, collapse=";"))
    mirs.weak <- data.frame(SYMBOL = rownames(mirs.weak), mir.weak = as.character(unlist(mirs.weak)),
                            stringsAsFactors=FALSE)
    mirs.weak.chr <- with(mi.gene[ind.weak,],
                          tapply(chr, Target.Gene, paste, collapse=";"))
    mirs.weak.chr <- data.frame(SYMBOL = rownames(mirs.weak.chr),
                                mir.weak.chr = as.character(unlist(mirs.weak.chr)))
    mirs.weak.pos <- with(mi.gene[ind.weak,],
                          tapply(start, Target.Gene, paste, collapse=";"))
    mirs.weak.pos <- data.frame(SYMBOL = rownames(mirs.weak.pos),
                                mir.weak.pos = as.character(unlist(mirs.weak.pos)))
    
    mirs.strong <- with(mi.gene[ind.strong,],
                      tapply(miRNA, Target.Gene, paste, collapse=";"))
    mirs.strong <- data.frame(SYMBOL = rownames(mirs.strong), mir.strong = as.character(unlist(mirs.strong)),
                            stringsAsFactors=FALSE)
    mirs.strong.chr <- with(mi.gene[ind.strong,],
                          tapply(chr, Target.Gene, paste, collapse=";"))
    mirs.strong.chr <- data.frame(SYMBOL = rownames(mirs.strong.chr),
                                mir.strong.chr = as.character(unlist(mirs.strong.chr)))
    mirs.strong.pos <- with(mi.gene[ind.strong,],
                          tapply(start, Target.Gene, paste, collapse=";"))
    mirs.strong.pos <- data.frame(SYMBOL = rownames(mirs.strong.pos),
                                mir.strong.pos = as.character(unlist(mirs.strong.pos)))
    
    genes <- merge(genes, mirs.weak, by="SYMBOL", all.x = TRUE, all.y=FALSE)
    genes <- merge(genes, mirs.weak.chr, by="SYMBOL", all.x = TRUE, all.y=FALSE)
    genes <- merge(genes, mirs.weak.pos, by="SYMBOL", all.x = TRUE, all.y=FALSE)
    genes <- merge(genes, mirs.strong, by="SYMBOL", all.x = TRUE, all.y=FALSE)
    genes <- merge(genes, mirs.strong.chr, by="SYMBOL", all.x = TRUE, all.y=FALSE)
    genes <- merge(genes, mirs.strong.pos, by="SYMBOL", all.x = TRUE, all.y=FALSE)
    cols <- c("mir.weak.chr", "mir.weak.pos", "mir.strong.chr", "mir.strong.pos")
    for (col in cols){
        genes[,col] <- as.character(genes[,col])
    }     

    ## SEX.IND : IS GENE ON X/Y? SEX.MIR*:DOES GENE INTERACT WITH A MIRNA ON X/Y?
    genes$sex.ind <- genes$sex.mir.weak.ind <-  genes$sex.mir.strong.ind <- 0
    xy.ind <- which(genes$chr %in% c("X","Y"))
    mir.xy.weak.ind <- grep("X|Y", genes$mir.weak.chr)
    mir.xy.strong.ind <- grep("X|Y", genes$mir.strong.chr)
    
    genes$sex.ind[xy.ind] <- 1
    genes$sex.mir.weak.ind[mir.xy.weak.ind] <- 1
    genes$sex.mir.strong.ind[mir.xy.strong.ind] <- 1
    
    ## COLLAPSE GENES CHRS / POSITIONS TOGETHER FOR SAME GENE SYMBOL ##
    dupe.ind.all <- which(duplicated(genes$ENSEMBL))
    dupes <- unique(genes$ENSEMBL[dupe.ind.all])
    col.names.local <- c("SYMBOL", "chr", "start","end", "intx",
                         "sex.ind", "sex.mir.weak.ind", "sex.mir.strong.ind")
    col.names.global <- c("mir.weak", "mir.weak.chr", "mir.weak.pos",
                          "mir.strong", "mir.strong.chr", "mir.strong.pos")

    for (dupe in dupes){
        dupe.ind <- which(genes$ENSEMBL == dupe)
        ## ARE ALL THE GENE SYMBOLS THE SAME? ##
        genes.differ <- length(unique(genes$SYMBOL[dupe.ind])) > 1
        for (col.name in col.names.local){
            if (col.name == "SYMBOL" & genes.differ==F) { next }
            col.ind <- which(names(genes) == col.name)            
            tmp <- genes[dupe.ind, col.ind]
            tmp <- paste(tmp, collapse=";")
            genes[dupe.ind, col.ind] <- tmp
        }

        if (genes.differ == F){ next }
        for (col.name in col.names.global){
            col.ind <- which(names(genes) == col.name)            
            tmp <- genes[dupe.ind, col.ind]
            tmp <- paste(as.character(tmp), collapse=")(")
            tmp <- paste0("(",tmp, ")")
            genes[dupe.ind, col.ind] <- tmp
        }        
    }
    
    genes <- genes[-dupe.ind.all,]

    return(genes)
}

get.gene.annotation <- function(rerun=F){
    library(biomaRt)
    gene.annot.file <- "/group/stranger-lab/askol/TCGA/Code/gene.annot.txt"

    if (rerun==F & file.exists(gene.annot.file)){

        print("Reading from file. Use rerun=T to regenerate file.")
        gene.info <- read.table(file = gene.annot.file, as.is=T, header=T)

    }else{

        ## GET ALL POSSIBLE GENE IDS IN TCGA DATA ##
        geneid <- system("awk '{print $1}' /group/stranger-lab/askol/TCGA/TCGA_Expression/TCGA-ACC.DATA.txt", intern=T)
        geneid <- gsub("\\.[0-9].*$","", geneid[-1])
        
        ## GET GENE ANNOTATION ##
        mart <- useMart("ensembl", host = "www.ensembl.org")
        mart <- useDataset(dataset = 'hsapiens_gene_ensembl',
                           mart = mart)        
        attribs = c("ensembl_gene_id", "hgnc_symbol",
            "chromosome_name", "start_position","end_position", "gene_biotype")
        gene.info <- getBM(attributes = attribs, filters = 'ensembl_gene_id',
                           values = geneid, mart = mart)
        names(gene.info) <- c("ENSEMBL", "SYMBOL", "chr", "start", "end",
                              "gene_biotype")
        
        ## REMOVE GENES WITH NO SYMBOL ##
        rm.ind <- which(is.na(gene.info$SYMBOL) | gene.info$SYMBOL == "")
        if (length(rm.ind)>0){
            gene.info <- gene.info[-rm.ind, ]
        }
        
        ## REMOVE DUPLICATES
        rm.ind <- which(duplicated(
            paste(gene.info$ENSEMBL, gene.info$chr, gene.info$SYMBOL,
                  sep=".")))
        if (length(rm.ind) > 0 ){
            gene.info <- gene.info[-rm.ind,]
        }
        
        write.table(file = gene.annot.file, gene.info, quote=F,
                    row.names=F, col.names=T)
    }    
    return(gene.info)
}

    

## ------ MAIN --------##
args <- commandArgs(TRUE)
project = args[1]
Create_Production_Data(project)
print(date())

finish.file = paste0(project,".finished")
system(paste0("touch ",finish.file))
