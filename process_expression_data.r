source("/group/stranger-lab/askol/TCGA/Code/miRNA-mRNA-Correlation_count_funcs.r")

DE_by_covariate <- function(project, cov="gender"){

    setwd("/group/stranger-lab/askol/TCGA/Expression_Analysis/")
    print(date())
    
    library("GenomicFeatures")
    library("RColorBrewer")
    library("biomaRt")
    library("AnnotationDbi")
    library("topGO")
    library("GO.db")
    library("org.Hs.eg.db")
    library(Homo.sapiens)
    library(limma)
    library(Glimma)
    library(edgeR)
    
    # ensembl = useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl"))

    EXPDIR = "/group/stranger-lab/askol/TCGA/TCGA_Expression/Limma_Voom_Use/"
    METADIR = "/group/stranger-lab/askol/TCGA/TCGA_MetaData/"
    ANALDIR = "/group/stranger-lab/askol/TCGA/Expression_Analysis/"
    GTEXDIR = "/group/gtex-group/DifferentialExpression/"
    gtexqFile =  paste0(GTEXDIR, "gender_de_limma_qval.tsv")
    gtexlfFile = paste0(GTEXDIR, "gender_de_limma_log2fc.tsv")

    ## GET DATA READY FOR LIMMA/VOOM ##
    count.auto <- readRDS(paste0(EXPDIR, project,"_mrna_autosome_DGEList.RDS"))
    count.sex <- readRDS(paste0(EXPDIR, project,"_mrna_sex_DGEList.RDS"))

    ## TEST AUTOSOMAL SNPS ##
    
    ## DETERMINE NUMBER PEERS TO BE CONDITIONED ON ##
    peer.names <- names(count.auto$samples)[grep("PEER",names(count.auto$samples))]
    no.sign.p <- -1
    no.sign <- 0
    vfit.p <- vfit <- c()
    efit.p <- efit <- c()
    peer.cnt <- 0
    
    ## FIND THE SET OF PEERS TO USE ##
    while(no.sign.p < no.sign & peer.cnt <= length(peer.names)){
        
        print(paste0("Usings peers: ", paste(peer.names[1:peer.cnt], collapse=",")))

        peer.cnt <- peer.cnt+1
        no.sign.p <- no.sign
        vfit.p <- vfit
        efit.p <- efit
        
        cmd <- paste("model.matrix(~0+group+",paste(peer.names[1:peer.cnt],collapse="+"),
                     ",count.auto$samples)")
        design <- eval(parse(text = cmd))
        colnames(design) <- gsub("group", "", colnames(design))       
        contr.matrix <- makeContrasts(Sex = male-female, levels = design)
        
        ## REMOVE HETEROSCEDASCITY ##
        ## In limma, linear modelling is carried out on the log-CPM values ##
        ## which are assumed to be normally distributed and the mean-variance ##
        ## relationship is accommodated using precision weights calculated by ##
        ## the voom function. ##
        v <- voom(count.auto, design, plot=FALSE)
        vfit <- lmFit(v, design)
        vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
        efit <- eBayes(vfit)
        
        no.sign <- sum(summary(decideTests(efit))[-2])
    }


    if (no.sign < no.sign.p){
        vfit <- vfit.p
        efit <- efit.p        
    }

    ## UPDATE PEER NAMES FOR GENES ON SEX CHROMOSOMES ##
    peer.names <- colnames(vfit$design)[grep("PEER",colnames(vfit$design))]
    cmd <- paste("model.matrix(~0+group+",paste(peer.names, collapse="+"),
                 ",count.sex$samples)")
    design <- eval(parse(text = cmd))
    colnames(design) <- gsub("group", "", colnames(design))       
    contr.matrix <- makeContrasts(Sex = male-female, levels = design)
    
    v.sex <- voom(count.sex, design, plot=FALSE)
    vfit.sex <- lmFit(v.sex, design)
    vfit.sex <- contrasts.fit(vfit.sex, contrasts=contr.matrix)
    efit.sex <- eBayes(vfit.sex)
    
    save(list = ls(all.names = TRUE), file=paste0(project,"_DE_limma.Rdata"))
      
}


## ------ MAIN --------##
args <- commandArgs(TRUE)
project = args[1]
DE_by_covariate(project)
print(date())

finish.file = paste0(project,".finished")
system(paste0("touch ",finish.file))

########################

