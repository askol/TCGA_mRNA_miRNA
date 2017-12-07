library("biomaRt")
library(UpSetR)
library(ggplot2)
library(reshape2)

source("/group/stranger-lab/askol/TCGA/Code/explore_expression_results_func.r")

GTEXDIR = "/group/gtex-group/DifferentialExpression/"
gtexqFile =  paste0(GTEXDIR, "gender_de_limma_qval.tsv")
gtexlfFile = paste0(GTEXDIR, "gender_de_limma_log2fc.tsv")
gtexpFile <- paste0(GTEXDIR, "gender_de_limma_pval.tsv")

setwd("/group/stranger-lab/askol/TCGA/Expression_Analysis/")

projects = read.table(file = "../TCGA_Target_projects.txt", header=TRUE, as.is=TRUE, sep="\t")
projects = projects[grep("TCGA",projects[,1]),1]



## get annotation for later ##
load("TCGA-ACC_DE.Rdata")
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

annot = result2[,c("ensembl","hgnc_symbol","chromosome_name","start_position",
    "end_position")]
annot <- as.data.frame(annot)

## GET GTEX DATA ##  
gtex.use = get_gtex(gtexpFile, gtexqFile, gtexlfFile)
gtex.use <- update.gtex(gtex.use)

## MIR GENES ##
mi.gene <- read.table("/group/stranger-lab/askol/TCGA/hsa_MTI-4.txt",
                      as.is=T, header=T, sep="\t")

mir.genes <- data.frame(mir.gene = unique(mi.gene$Target.Gene[-grep("Weak", mi.gene$Support.Type)]))
mir.genes <- merge(mir.genes, resultTable[,c("hgnc_symbol", "ensembl_gene_id")], by=1,
                   keep.x=T, keep.y=F)

## SINGLE SEX CANCERS ##
skip.cancers <- c("TCGA-BRCA","TCGA-OV", "TCGA-CESC", "TCGA-PRAD", "TCGA-TGCT", "TCGA-UCS", "TCGA-UCEC")
## skip.tmp <- "TCGA-LUAD"

sex.de.tbl = c()
de.genes = list()

## STORE RESULTS VARIABLE ##
## store <- list()

for (project in projects){
    
    print(paste0("Working on project ",project))
    if (project %in% c(skip.cancers)){
        print("Skipping single sex cancer. . .")
        next
    }
    
    ## file = paste0(project, "_DE.Rdata")
    file = paste0(project, "_result2.RData")
    if (!file.exists(file)){ next }

    load(file)
  
    ## results2 is the relavent data.frame ##
    ## GET GENES WITH QVALUES OF 0.05 OR LESS ##
    result2$sig = 1*(result2$padj <= 0.05)
    result2 = as.data.frame(result2)
    result2 = merge(result2, gtex.use, by = "ensembl", all.x=T, all.y=F)
    result2 = merge(result2, mir.genes, by.x="ensembl", by.y="ensembl_gene_id", all.x=T, all.y=F)
    result2$mir.gene.ind <- (is.na(result2$mir.gene)==FALSE)*1

    ## FILL SIGN.DE GENES FOR GRAPHING OVERLAP ##
    #if (length(sign.DE.genes) == 0){
    #    sign.DE.genes <- result2[,c("ensembl", "sig")]
    #    names(sign.DE.genes)[2] = paste0(project,".sig")
    #}else{
    #    sign.DE.genes <- merge(sign.DE.genes, result2[,c("ensembl","sig")], by=1, all=T)
    #    names(sign.DE.genes)[ncol(sign.DE.genes)] = paste0(project,".sig")
    #}
            
    ## LOOK UP MATCHING NORMAL TISSUE(S)
    top.tiss <- cancer_tissue_lookup(project)
    
    col.ind = intersect(grep(paste(top.tiss, collapse="|"), names(result2)),
        grep("qval",names(result2)))
    result2$sig.norm = 1*(rowSums(as.matrix(result2[,col.ind])<=.3, na.rm=T) > 0)
    all.na.ind <- which(rowSums(is.na(as.matrix(result2[,col.ind]))) == length(col.ind))
    result2$sig.norm[all.na.ind] = NA

    ## what proportion of genes are DE in cancer wrt sex?
    prop.de.cancer <- mean(result2$sig, na.rm=T)
    no.de.cancer <- sum(result2$sig, na.rm=T)
    no.de.normal <- sum(result2$sig.norm, na.rm=T)
    
    ## what proportion of de cancer genes found in sim. tissue ##
    no.de.norm.canc <- sum(result2$sig.norm == 1 & result2$sig == 1, na.rm=T)
    no.de.canc.not.norm <- no.de.cancer - no.de.norm.canc
    no.de.cancer.adj.norm = sum(result2$sig[ -which(result2$sig.norm==1)], na.rm=T)

    prop.red.de.cancer <- 1 - (no.de.cancer.adj.norm)/(no.de.cancer)
    prop.de.cancer.adj.norm = mean(result2$sig[ -which(result2$sig.norm==1)],
        na.rm=T)
    
    prop.uniq.to.cancer = 1 - prop.red.de.cancer

    ## what proportion of de cancer genes are in sim.tissue and mir target ##
    no.de.cancer.mir <- sum(result2$mir.gene.ind & result2$sig, na.rm=T)
    no.de.normal.cancer.mir <- sum(result2$mir.gene.ind & result2$sig & result2$sig.norm, na.rm=T)
    no.de.cancer.not.norm.mir <- no.de.cancer.mir - no.de.normal.cancer.mir 
    
    ## STORE FOUND GENES ##
    ind = which(result2$sig == 1 & result2$sig.norm==0)
    genes.tmp <-  unique(result2$ensembl[ind])
    out.tmp <- merge(result2[ind,c("ensembl","log2FoldChange","pvalue","padj")],
                     annot, by="ensembl", all.x=T, all.y=F)
    tmp <- merge(result2[,c("ensembl","sig","sig.norm")],
                 annot[,c("ensembl","hgnc_symbol")], by="ensembl",
                 all.x=T, all.y=F)
    tmp$sig[result2$sig.norm==1] = NA
    names(tmp)[2] = project
    de.genes[[project]] <- tmp[,c("ensembl",project)]
    
    file = paste0(project,"_DM_results.txt")
    write.table(file = file, out.tmp, quote=F, col.names=T, row.names=F)
    
    sex.de.tbl = rbind(sex.de.tbl,
        c(project, no.de.cancer, prop.de.cancer, no.de.normal,
          prop.de.cancer.adj.norm, no.de.norm.canc, no.de.canc.not.norm,
          prop.red.de.cancer, prop.uniq.to.cancer,
          no.de.cancer.mir, no.de.normal.cancer.mir, no.de.cancer.not.norm.mir,
          paste(top.tiss,collapse="|")))

    if (0){
        ## CREATE QQ PLOTS OF FOR NORMAL DE GENES IN CANCER AND VISE VERSA ##
        gtex.ind <- grep(paste(gsub("\\.+qval","",top.tiss), collapse="|"), names(gtex.use))
        plot.p.dist.canc.v.norm(project, gtex.use[,c(1, gtex.ind)], result2[,c("ensembl","pvalue","padj")]) 
        
        ## store[[project]] <- result2
    }
    
    sizes <- sapply(ls(), function(n) object.size(get(n)),
                    simplify = FALSE);
    print(sapply(sizes[order(as.integer(sizes))],
                 function(s) format(s, unit = 'auto')))
}

colnames(sex.de.tbl) = c("project", "No.DE.canc.genes",
            "Prop.DE.canc.genes", "No.DE.in.norm",
            "prop.DE.cancer.adj.norm",
            "No.DE.Canc.and.Norm", "No.DE.Canc.not.norm",
            "Prop.Red.DE.Cancer","Prop.uniq.Canc",
            "No.DE.cancer.mir", "No.DE.normal.cancer.mir",
            "No.DE.cancer.not.normal.mir","Tissue.Norm")

write.table(file = "Sex.DE.table.txt", sex.de.tbl, row.names=F,
            col.names=T, quote=F)

time.stamp <- strftime(Sys.time(), "%Y%m%d%H%M%S")
save.image(paste0("explore_expression_results_",time.stamp,".RData"))
q()

## R CMD BATCH explore_expression_results.r explore_expression_results.out&
sex.de.tbl <- as.data.frame(sex.de.tbl)
sex.de.tbl[, 2:7] <- apply(sex.de.tbl[, 2:7], 2, function(x) as.numeric(as.character(x)) ) 

mean(sex.de.tbl$Prop.red.DE.cancer.genes[-c(1:2)])
range(sex.de.tbl$Prop.red.DE.cancer.genes[-c(1:2)])

mean(sex.de.tbl$prop.DE.cancer.adj.norm[-c(1:2)])
range(sex.de.tbl$prop.DE.cancer.adj.norm[-c(1:2)])

mean(sex.de.tbl$Prop.uniq.to.cancer[-c(1:2)])
range(sex.de.tbl$Prop.uniq.to.cancer[-c(1:2)])


### PLOTS ##
## Determine distribution genes shared across cancers ##
genes <- de.genes[[1]]
genes <- genes[which(genes[,2] == 1),]
for (project in names(de.genes)[-1]){
    print(project)
    tmp <- de.genes[[project]]
    tmp <- tmp[which(tmp[,2]==1),]
    genes <- merge(genes, tmp, by = "ensembl", all=T)
}
genes[is.na(genes)] = 0
plot.de.gene.overlap(genes)

## PLOT BAR PLOT SHOWING NUMBER OF DE GENES IN EACH CANCER STACKED FOR SHARING
## IN OTHER CANCER ##
genes.ol <- genes[genes$ensembl %in% mir.genes[,2],] 
plot.no.de.mir.gene.ol(genes.ol)

## PLOT DISTRIBUTION OF THE NUMBER OF DE GENES PER CANCER
## -- INCLUDE NUMBER DISTINCT TO CANCER VS NORMAL
plot.no.de.genes(sex.de.tbl)

## PLOT DISTN OF NO DE GENES PER CANCER FOR ONLY GENES TARGETTING MIRS
plot.no.de.mir.genes(sex.de.tbl)


q()
##### --- functions

all.genes = unique(unlist(de.genes))
canc.genes = matrix(0, length(all.genes), length(de.genes)-1)
canc.genes[,1] = all.genes
for (i in 3:length(de.genes)){

    col = i-1

    ind = match(de.genes[[i]], all.genes)
    canc.genes[ind, col] = 1

}
canc.genes = as.data.frame(canc.genes)
canc.genes[,-1] = apply(canc.genes[,-1], 2, function(x) as.numeric(as.character(x)))
table(rowSums(canc.genes[,-1]))
    
