


get_gtex <- function(pFile, qFile, lfFile){
    
    data = read.table(qFile, as.is=TRUE, header=T, sep="\t")
    names(data)[1] = "ensembl"
    names(data)[-1] = paste(names(data)[-1], "qval", sep=".")
    data$ensembl = gsub("\\.\\d+", "", data$ensembl)
    
    lf = read.table(lfFile, as.is=TRUE, header=T, sep="\t")
    names(lf)[1] = "ensembl"
    names(lf)[-1] = paste(names(lf)[-1], "lf", sep=".")
    lf$ensembl = gsub("\\.\\d+", "", lf$ensembl)

    data = merge(data, lf, by=1, all.x=TRUE)
    
    pval = read.table(pFile, as.is=TRUE, header=T, sep="\t")
    names(pval)[1] = "ensembl"
    names(pval)[-1] = paste(names(pval)[-1], "pval", sep=".")
    pval$ensembl = gsub("\\.\\d+", "", pval$ensembl)
    
    data = merge(data, pval, by=1, all.x=TRUE)
    
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


update.gtex <- function(gtex){
    ind = grep("qval", names(gtex))
    nms = names(gtex)[ind]
    for (nm in nms){
        nw.nm = gsub("qval","q.lt.3", nm)
        eval(parse(text = paste0("gtex$",nw.nm,"=1*(gtex$",nm,"<.3)")))
    }
    return(gtex)
}


calc.concord <- function(g, r, qthresh = .3){

    inds = grep("q.lt.3", names(g))
    r = as.data.frame(r)
    out = c()
    
    for (tis.ind in inds){

        nm = names(g)[tis.ind]
        if (length(grep("Breast",nm))){ next }
        tis = gsub(".q.lt\\.3", "", perl=T, nm)

        r$q.lt = 1*(r$padj <= qthresh)
        
        tmp = merge(g[,c("ensembl",nm)], r[,c("ensembl","q.lt")],
            by = "ensembl")

        tbl = table(tmp[,2:3])
        or = tbl[1,1]*tbl[2,2]/tbl[1,2]/tbl[2,1]
        p = chisq.test(tbl)$p.value

        out = rbind(out, c(tis, or, p))

    }

    out = data.frame(out)
    colnames(out) = c("tissue","OR","P")
    out[,-1] = apply(out[,-1], 2, function(x) as.numeric(as.character(x)))
    out[,1] = as.character(out[,1])
    return(out)

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
                     "TCGA-KIRC" = c("Liver","Small.Intestine...Terminal.Ileum"),
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

     
plot.no.de.genes <- function(count.tbl){

    plot.name <- "No_DE_genes_by_project.pdf"
    ## shape count table ##
    dat <- as.data.frame(count.tbl)
    dat <- count.tbl[,c("project","No.DE.Canc.and.Norm","No.DE.Canc.not.norm")]
    dat[,-1] <- apply(dat[,-1],2, function(x) as.numeric(as.character(x)))
    
    ord <- order(rowSums(dat[,-1]), decreasing = TRUE)
    dat <- dat[ord,]
                 
    dat <- melt(dat, id.vars="project", variable.name = "gene.type", value.name="no.genes")    

    ## GET NUMBER OF SAMPLES PER PROJECT ##
    ## ADJUST PROJECT NAMES TO INCLUDE SAMPLE SIZE
    header <- read.table(file = "/group/stranger-lab/askol/TCGA/TCGA_MetaData/Sex_table.txt",as.is=T, nrow=1)
    samp.count <- read.table(file = "/group/stranger-lab/askol/TCGA/TCGA_MetaData/Sex_table.txt",as.is=T, skip=1)
    names(samp.count) <- c("project",header)

    samp.count$m.count <- rowSums(samp.count[,c("m.male","m.female")])
    samp.count <- samp.count[samp.count$project %in% dat$project,]
    samp.count <- samp.count[match(dat$project, samp.count$project), ]
    samp.count$label <- paste0(samp.count$project, " (N = ", samp.count$m.count,")")
    dat$project <- samp.count$label
    dat$project = factor(dat$project, levels = dat$project[1:(nrow(dat)/2)])
    pdf(plot.name, width = 8, height=6)
    p <- ggplot(dat) +
        geom_col(aes(x = project, y = no.genes, fill = gene.type), position = "stack") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))

    print(p)
    dev.off()
}

## SAME PLOT AS ABOVE BUT FOR GENES TARGETED BY MIR ##
plot.no.de.mir.genes <- function(count.tbl){

    plot.name <- "No_DE_mir.genes_by_project.pdf"
    ## shape count table ##
    dat <- as.data.frame(count.tbl)
    dat <- count.tbl[,c("project","No.DE.normal.cancer.mir", "No.DE.cancer.not.normal.mir")]
    dat <- as.data.frame(dat)
    dat[,-1] <- apply(dat[,-1],2, function(x) as.numeric(as.character(x)))
    
    ord <- order(rowSums(dat[,-1]), decreasing = TRUE)
    dat <- dat[ord,]
                 
    dat <- melt(dat, id.vars="project", variable.name = "gene.type", value.name="no.genes")    

    ## GET NUMBER OF SAMPLES PER PROJECT ##
    ## ADJUST PROJECT NAMES TO INCLUDE SAMPLE SIZE
    header <- read.table(file = "/group/stranger-lab/askol/TCGA/TCGA_MetaData/Sex_table.txt",as.is=T, nrow=1)
    samp.count <- read.table(file = "/group/stranger-lab/askol/TCGA/TCGA_MetaData/Sex_table.txt",as.is=T, skip=1)
    names(samp.count) <- c("project",header)

    samp.count$m.count <- rowSums(samp.count[,c("m.male","m.female")])
    samp.count <- samp.count[samp.count$project %in% dat$project,]
    samp.count <- samp.count[match(dat$project, samp.count$project), ]
    samp.count$label <- paste0(samp.count$project, " (N = ", samp.count$m.count,")")
    ## REMOVE TCGA ##
    samp.count$label <- gsub("TCGA-","", samp.count$label)
    
    dat$project <- samp.count$label
    dat$project = factor(dat$project, levels = dat$project[1:(nrow(dat)/2)])
    pdf(plot.name, width = 8, height=6)
    p <- ggplot(dat) +
        geom_col(aes(x = project, y = no.genes, fill = gene.type), position = "stack") +
            theme_minimal() +
                theme(text=element_text(size=14),
                      axis.text.x = element_text(angle = 90, hjust = 1))+
                          ylab("Number DE Genes")
                           
    print(p)
    dev.off()
}

## 
plot.no.de.mir.gene.ol <- function(genes,
                                   plot.name = "No_DE_mir.genes_by_project_ol.pdf")
    {

        genes <- genes[,-1]
        tbl <- c()
        for (col in 1:ncol(genes)){
            
            no.genes <- sum(genes[,col], na.rm=T)
            no.genes.ol <- sum(genes[,col]==1 & rowSums(genes[,-col],
                                        na.rm=T)>0)
            tbl <- rbind(tbl, c(names(genes)[col], "unique", no.genes-no.genes.ol))
            tbl <- rbind(tbl, c(names(genes)[col], "shared", no.genes.ol))
        }
        tbl <- as.data.frame(tbl)
        names(tbl) <- c("project", "gene.type", "no.genes")
        tbl$no.genes <- as.numeric(as.character(tbl$no.genes))   
        
        header <- read.table(file = "/group/stranger-lab/askol/TCGA/TCGA_MetaData/Sex_table.txt",as.is=T, nrow=1)
        samp.count <- read.table(file = "/group/stranger-lab/askol/TCGA/TCGA_MetaData/Sex_table.txt",as.is=T, skip=1)
        names(samp.count) <- c("project",header)
        
        samp.count$m.count <- rowSums(samp.count[,c("m.male","m.female")])
        samp.count <- samp.count[samp.count$project %in% tbl$project,]
        samp.count <- samp.count[match(tbl$project, samp.count$project), ]
        samp.count$label <- paste0(samp.count$project, " (N = ", samp.count$m.count,")")
        ## REMOVE TCGA ##
        samp.count$label <- gsub("TCGA-","", samp.count$label)
        
        tbl$project <- samp.count$label
        tmp <-dcast(tbl, project~gene.type)
        levels <- tmp$project[order(rowSums(tmp[,-1]), decreasing=TRUE)]
        tbl$project = factor(tbl$project, levels =levels)
        
        pdf(plot.name, width = 8, height=6)
        p <- ggplot(tbl) +
            geom_col(aes(x = project, y = no.genes, fill = gene.type), position = "stack") +
                theme_minimal() +
                    theme(text=element_text(size=14),
                      axis.text.x = element_text(angle = 90, hjust = 1))+
                          ylab("Number DE Genes")
        
        print(p)
        dev.off()
        
    }

plot.p.dist.canc.v.norm <- function(project, norm.ps, canc.ps){

    plot.file <- paste0(project,"qq_plots.pdf")
    
    ## DETERMINE NORMAL TISSUES ##
    tissues <- names(norm.ps)[grep("qval",names(norm.ps))]
    sign.canc <- canc.ps$ensembl[canc.ps$padj < 0.05]

    pdf(plot.file)
    
    for (tissue in tissues){
        
        tissue <- gsub("\\.qval","",tissue)
        tissue.ind <- grep(tissue, names(norm.ps))
        tissue.p.ind <- intersect(tissue.ind, grep("pval", names(norm.ps)))
        tissue.q.ind <- intersect(tissue.ind, grep("qval", names(norm.ps)))
        sign.norm <- norm.ps$ensembl[norm.ps[,tissue.q.ind] <= .05]
        tissue.txt <- gsub("\\.\\.\\..*", "",tissue)
        title.pre <- paste0(project," : ",tissue.txt)
        
        ## QQ PLOT OF ALL P-VALUES FROM CANCER DE ANALYSIS ##
        title <- paste(title.pre, "- All genes")
        p <- qq(canc.ps$pvalue, title=title)
        print(p)
        
        ## QQ PLOT OF P-VALUES OF CANC TISSUE WITH NORMAL DE GENES REMOVED ##
        title <- paste(title.pre, "- Cancer|GTEx.DE.removed")
        p.cancer.g.norm <- canc.ps$pvalue[!canc.ps$ensembl %in% sign.norm]
        p<-qq(p.cancer.g.norm, title=title)
        print(p)
        
        ## QQ PLOT OF P-VALUES OF CANC TISSUE FOR NORMAL NON-DE GENES ##
        title <- paste(title.pre, "- Cancer|GTEx.DE")
        p.cancer.g.norm <- canc.ps$pvalue[canc.ps$ensembl %in% sign.norm]
        p<-qq(p.cancer.g.norm, title=title)
        print(p)
        
        ## QQ PLOT OF P-VALUES OF NORM TISSUE FOR CANCER DE GENES ##
        title <- paste(title.pre, "- GTEx|cancer.DE")
        p.norm.g.cancer <- norm.ps[norm.ps$ensembl %in% sign.canc, tissue.p.ind]
        p <- qq(p.norm.g.cancer, title=title)
        print(p)
    
    }

    dev.off()

    print(paste("Created plot",plot.file))

}
    
    
plot.de.gene.overlap <- function(genes){

    ## upsetr is expecting project as columns and genes as rows ##
    #genes <- t(genes)

    ## GET NUMBER OF SAMPLES PER PROJECT ##
      ## GET NUMBER OF SAMPLES PER PROJECT ##
    header <- read.table(file = "/group/stranger-lab/askol/TCGA/TCGA_MetaData/Sex_table.txt",as.is=T, nrow=1)
    samp.count <- read.table(file = "/group/stranger-lab/askol/TCGA/TCGA_MetaData/Sex_table.txt",as.is=T, skip=1)
    names(samp.count) <- c("project",header)

    samp.count$m.count <- rowSums(samp.count[,c("m.male","m.female")])
    samp.count <- samp.count[samp.count$project %in% names(genes),]
    samp.count <- samp.count[match(names(genes)[-1], samp.count$project), ]
    samp.count$label <- paste0(samp.count$project, " (N = ", samp.count$m.count,")")
    names(genes)[-1] <- samp.count$label
    
    pdf(file = "Intersection_DM_genes.pdf", width=10, height=7)
    upset(genes, sets = names(genes)[-1], sets.bar.color = "#56B4E9",
          order.by = "freq")
    dev.off()

    print("Created plot Intersection_gene.pdf")
}

qq <- function(pvector, title="Quantile-quantile plot of p-values", spartan=F) {
    # Thanks to Daniel Shriner at NHGRI for providing this code for creating expected and observed values
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    # you could use base graphics
    #plot(e,o,pch=19,cex=0.25, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)),ylim=c(0,max(e)))
    #lines(e,e,col="red")
    #You'll need ggplot2 installed to do the rest
    plot=qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o))) + geom_abline(intercept=0,slope=1, col="red")
    plot=plot+ggtitle(title)
    plot=plot+scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))
    plot=plot+scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))
    if (spartan) plot=plot+opts(panel.background=theme_rect(col="grey50"), panel.grid.minor=theme_blank())
    return(plot)
}
