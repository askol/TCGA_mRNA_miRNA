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
        vfit.p <- vfit
        efit.p <- efit
        
        cmd <- paste("model.matrix(~0+group+",paste(peer.names[1:peer.cnt],collapse="+"),
                     ",count.auto$samples)")
        design <- eval(parse(text = cmd))
        colnames(design) <- gsub("group", "", colnames(designs[[i]]))       
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


    if (no.sign.p < no.sign){
        vfit <- vfit.p
        efit <- efit.p        
    }

    ## UPDATE PEER NAMES FOR GENES ON SEX CHROMOSOMES ##
    peer.names <- colnames(vfit$design)[grep("PEER",colnames(vfit$design))]
        
    

    print("Creating DataSet. . .")
    dataset = DESeqDataSetFromMatrix(countData = countData, colData = meta.data.use,
        design = ~gender)

    print("Running DESeq. . .")
    
    dds <- DESeq(dataset)
    GeneNamesAll = names(dds)
    nrow(dds)
    ## REMOVE GENES THAT HAVE VERY LOW COUNTS ##
    dds <- dds[ rowSums(counts(dds)) > ncol(dds)/2, ]
    nrow(dds)

    ## TRY TO IDENTIFY TECHNICAL ARTIFACTS ##
    print("Running svaseq")
    dat  <- counts(dds, normalized = TRUE)
    idx  <- rowMeans(dat) > 1
    dat  <- dat[idx, ]
    mod  <- model.matrix(~ gender, colData(dds))
    mod0 <- model.matrix(~   1, colData(dds))
    svseq <- svaseq(dat, mod, mod0, n.sv = 2)

    ddssva <- dds
    ddssva$SV1 <- svseq$sv[,1]
    ddssva$SV2 <- svseq$sv[,2]
    design(ddssva) <- ~ SV1 + SV2 + gender

    print("Running DESeq ond ddssva")
    dds2 = DESeq(ddssva)
    GeneNamesAll = names(dds2)
    
    print("Running results. . .")
    result2 <- results(dds2, contrast=c(cov,levels(meta.data.use[[cov]])))
    
    rownames(result2) =  gsub("\\.\\d+$", "", rownames(result2))
    result2$ensembl = rownames(result2)

    print("Writing file and saving R data")
    write.table(file = paste0(project,"_DE_results.txt"), result2, quote=F,
                row.names=TRUE, col.names=TRUE)
    save(list = ls(all.names = TRUE), file=paste0(project,"_DE.Rdata"))
    
    ## Annotate results significant results##
    gene_set = result2$ensembl
    mart <- useDataset(dataset = "hsapiens_gene_ensembl",
                       mart = useMart("ensembl",
                           host    = "www.ensembl.org"))

    attribs = c("ensembl_gene_id", "entrezgene", "hgnc_symbol",
        "chromosome_name", "start_position","end_position")
    resultTable <- getBM(attributes = attribs, filters = "ensembl_gene_id",
                         values = gene_set, mart = mart)

    mtch = match(resultTable$ensembl_gene_id, result2$ensembl)
    
    for (attrib in attribs){

        eval(parse(text = paste0("result2$", attrib, "=NA")))
        eval(parse(text = paste0("result2$", attrib, "[mtch] = resultTable$",attrib)))
    }

    resultSig = result2[which(result2$padj <= 0.05), ]    
    
    AllGenes = result2$pvalue
    names(AllGenes) = result2$entrezgene
    ## REMOVE GENES FOUND DIFFERENTIALLY EXPRESSED IN GTEX BLOOD OR LCLS ##
    AllGenes = AllGenes[!names(AllGenes) %in% gtexSign]
    InterestingGenes = function(AllGenes) { return(AllGenes <= 0.05) }
   
    GSEA = run.gsea(AllGenes, InterestingGenes)

    print("Writing R data again")
####
    save(list = ls(all.names = TRUE), file=paste0(project,"_DE.Rdata"))
      
}

run.gsea <- function(AllGenes, InterestingGenes){

    GSEA.rslts = list()
    
    goLevels = c("BP", "MF", "CC")

    for (i in 1:3){

        goLevel = goLevels[i]
        GOdata <- new("topGOdata", ontology = goLevel, allGenes = AllGenes,
                      geneSel = InterestingGenes,
                      nodeSize = 5, description = paste0(project,"-",goLevel),
                      annot = annFUN.org,
                      mapping="org.Hs.eg.db",
                      ID="entrez" )
        
        resultFisher <- runTest(GOdata, algorithm="classic", statistic = "fisher")
        resultKS <- runTest(GOdata, algorithm="classic", statistic = "ks")
        resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
        

        GSEA.rslts[[i]] <- GenTable(GOdata, classicFisher = resultFisher,
                                classicKS = resultKS, elimKS = resultKS.elim,
                                orderBy = "elimKS", ranksOf = "classicFisher",
                                topNodes = 1000)

        out.file = paste0(project,"-GSEA-",goLevel,".out")
        write.table(file = out.file, GSEA.rslts[[i]], quote=FALSE, col.names=TRUE)
        print(paste("Wrote file", out.file))
    }
    return(GSEA.rslts)
}

annot.results <- function(res){

    gene_set = rownames(res)
    gene_set = gsub("\\.\\d+", "", gene_set)
    rownames(res) = gene_set
    
    mart <- useDataset(dataset = "hsapiens_gene_ensembl", 
                       mart = useMart("ensembl",
                           host = "www.ensembl.org"))

    attribs = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position")
    
    resultTable <- getBM(attributes = attribs,
                     filters    = "ensembl_gene_id",
                     values     = gene_set, 
                     mart       = mart)

    mtch = match(resultTable$ensembl_gene_id, rownames(res))
    for (attrib in attribs){
        txt = paste0("res$",attrib," = NA")
        eval(parse(text = txt))
        txt = paste0("res$",attrib,"[mtch] = resultTable$",attrib)
        eval(parse(text = txt))
    }

    ord = order(res$pvalue)
    res = res[ord,]
    
    return(res)
}

get_gtex <- function(qFile, lfFile){
    
    data = read.table(qFile, as.is=TRUE, header=T, sep="\t")
    names(data)[1] = "ensembl"
    names(data)[-1] = paste(names(data)[-1], "qval", sep=".")
    data$ensembl = gsub("\\.\\d+", "", data$ensembl)
    
    lf = read.table(lfFile, as.is=TRUE, header=T, sep="\t")
    names(lf)[1] = "ensembl"
    names(lf)[-1] = paste(names(lf)[-1], "lf", sep=".")
    lf$ensembl = gsub("\\.\\d+", "", lf$ensembl)
    
     data = merge(data, lf, by=1, all.x=TRUE)
    
    return(data)
}

get_gtex_sign <- function(gtex, tissues = c("EBV.transformed", "Whole.Blood"), pThresh = .05){
    
    ind = grep(paste(tissues, collapse="|"), names(gtex))
    ind = intersect(ind, grep("qval",names(gtex)))
    
    if (length(ind) != length(tissues)){
        print("One or more tissues not found")
        stop()
    }

    sign.genes.ind = which(rowSums( gtex[, ind] <= pThresh ) > 0)
    sign.genes = gtex[sign.genes.ind, 1]

    return(sign.genes)
}

cancer_tissue_lookup <- function(project){

    lu.table <- list("TCGA-LAML"=c("Cells...EBV.transformed.lymphocytes","Whole.Blood"),
                     "TCGA-ACC" = c("Adrenal.Gland"),
                     "TCGA-BLCA" = c("Stomach","Spleen"),
                     "TCGA-LGG" = c("Brain...Cortex"),
                     "TCGA-CHOL" = c("Stomach","Spleen"),
                     "TCGA-COAD" = "Colon...Sigmoid",
                     "TCGA-ESCA" = "Esophagus...Mucosa",
                     "TCGA-GBM" = c("Brain...Cortex","Brain...Spinal.cord..cervical.c.1"),
                     "TCGA-HNSC" = c("Minor.Salivary.Gland","Esophagus...Mucosa"),
                     "TCGA-KICH" = c("Liver","Small.Intestine...Terminal.Ileum"),
                     "TCGA-KIRK" = c("Liver","Small.Intestine...Terminal.Ileum"),
                     "TCGA-KIRP" = c("Liver","Small.Intestine...Terminal.Ileum"),
                     "TCGA-LIHC" = "Liver",
                     "TCGA-LUAD" = "Lung",
                     "TCGA-LUSC" = "Lung",
                     "TCGA-DLBC" = c("Cells...EBV.transformed.lymphocytes","Whole.Blood"),
                     "TCGA-MESO" = "Lung",
                     "TCGA-PAAD" = "Pancreas",
                     "TCGA-PCPG" = "Adrenal.Gland",
                     "TCGA-READ" = "Colon",
                     "TCGA-SARC" = c("Muscle...Skeletal","Adipose...Subcutaneous",
                     "mAdipose...Visceral", "Whole.Blood"),
                     "TCGA-SKCM" = "Skin...Not.Sun.Exposed..Suprapubic",
                     "TCGA-STAD" = "Stomach",
                     "TCGA-TGCT" = "Testis",
                     "TCGA-THYM" = c("Cells...EBV.transformed.lymphocytes","Whole.Blood"),
                     "TCGA-THCA" = "Thyroid",
                     "TCGA-UVM" = "Skin...Not.Sun.Exposed..Suprapubic"
                     )

    return(lu.table[[project]])
}

                     


## ------ MAIN --------##
args <- commandArgs(TRUE)
project = args[1]
DE_by_covariate(project)
print(date())

finish.file = paste0(project,".finished")
system(paste0("touch ",finish.file))

########################






     
#### OLD PLOTTING CODE ####
if (0){
    plot.name = paste0(ANALDIR, project,"-ma.pdf")
    print(paste0("Creating plot ",plot.name))
    pdf(file = plot.name)
    plotMA(result, main=paste0("Condition: ",lvls[1], " vs. ", lvls[2]), ylim=c(-5,5))
    dev.off()
    print(paste0("Wrote ",plot.name))

    print("Running rlogTransformation. . .")
    rld <- rlogTransformation(dds, blind=TRUE)

    ## SAVE DATA ##
    out.file = paste0(project,".RData")
    save.image(file = out.file)
    
    plot.name = paste0(ANALDIR, project, "-pca.pdf")
    print(paste0("Creating plot ", plot.name))
    pdf(file = plot.name)
    plotPCA(rld, intgroup = cov)
    dev.off()
    print(paste0("Wrote ",plot.name))
    
    ## Plot counts for 10 genes with the most significatn p-values
    ord = order(result$padj)
    plot.name = paste0(ANALDIR, project,"-top10.pdf")
    pdf(file = plot.name)
    print(paste0("Creating plot ",plot.name))
    for (ind in ord[1:10]){
        plotCounts(dds, gene=ind, intgroup=cov, pch = 19)
    }
    dev.off()
    print(paste0("Wrote ",plot.name))
    
    ## Display these top genes’ normalized counts in a heatmap, and cluster samples by similarity:
    n = 250 
    resOrdered <- result[order(result$padj),]
    topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,],
                        resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
    
    plot.name = paste0(ANALDIR, project, "-heatmap.pdf")
    print(paste0("Creating plot ", plot.name))
    pdf(file = plot.name)
    hmcol <- brewer.pal(11,'RdBu')
    nCounts <- counts(dds, normalized=TRUE)
    heatmap(as.matrix(nCounts[ row.names(topResults), ]), Rowv = NA, col = hmcol, mar = c(10,2))
    
    ## Examine sample clusters that arise from the top 25 and bottom 25 genes are used:
    m = 25
    heatmap(as.matrix(nCounts[ row.names(topResults)[c(1:m,(n-m+1):n)], ]), Rowv = NA, col = hmcol, mar = c(10,2))
    dev.off()

    print(paste0("Wrote ",plot.name))
    
    ##  Write results to file
    out.file = paste0(ANALDIR, project, "_results.txt")
    print(paste0("Writing to ",out.file))
    write.table(file = out.file, result, quote = FALSE, sep = '\t', col.names=T)
    print(paste0("Wrote ",out.file))

    finish.file = paste0(ANALDIR, project,".finished")
    system(paste0("touch ",finish.file))
}
