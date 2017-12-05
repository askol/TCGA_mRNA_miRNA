make.mi.map <- function(project){

    METADIR <- "/group/stranger-lab/askol/TCGA/TCGA_MetaData"
      
    meta.miRNA.file <- paste0(METADIR,"/",project,"_miRNA_meta.txt")
      
    mi.map <- read.table(file=meta.miRNA.file, header=T, as.is=T, sep="\t")
    names(mi.map) <- gsub("file_id", "file_id_mi", names(mi.map))
    names(mi.map) <- gsub("cases.0.samples.0.","", names(mi.map))
    
    return(mi.map)
  }

make.m.map <- function(project){

     METADIR <- "/group/stranger-lab/askol/TCGA/TCGA_MetaData"
     meta.mRNA.file <- paste0(METADIR,"/",project,"_mRNA_meta.txt")
     
     m.map <- read.table(file=meta.mRNA.file, header=T, as.is=T, sep="\t")
     
     names(m.map) <- gsub("file_name", "file_name_m", names(m.map))
     m.map$file_name_m = gsub("\\..*","",m.map$file_name_m)
     
     names(m.map) <- gsub("cases.0.samples.0.","", names(m.map))
     
     return(m.map)
 }

make.map.files <- function(project){
    
    m.map <- make.mi.map.file(project)
    mi.map <- make.m.map.file(project)
    
    return(list(m.map = m.map, mi.map = mi.map))

}

get.gene.info <- function(){

    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    gene.info <- getBM(attributes=c('ensembl_gene_id',
                           'ensembl_transcript_id','hgnc_symbol','chromosome_name',
                           'start_position','end_position', 'gene_biotype'), 
                       mart = ensembl)
    gene.info$id <- paste(gene.info$ensembl_gene_id, gene.info$hgnc_symbol,
                          sep=".")
    gene.info <- gene.info[!duplicated(gene.info$id),]
    return(gene.info)
}

if (0){

    dupe.ids <- unique(m.map$cases.0.samples.0.sample_id[duplicated( m.map$cases.0.samples.0.sample_id)])
    
    for (id in dupe.ids){
        print(m.map[m.map$cases.0.samples.0.sample_id == id,])
    }
}

get.data <- function(m.file, mi.file, m.map, mi.map, gene.info){

    ## GET MI_RNA DATA ##
    mi <- get.mi.data(mi.file, mi.map)
    
    ## GET MRNA DATA ##
    m <- get.m.data(m.file, m.map, gene.info)
    return(list(m = m, mi = mi))
}    
    
    
get.mi.data <- function(file, mi.map){

    head <- read.table(file = file, as.is=T, nrow=1, header=F)
    mi <- read.table(file = file, as.is=T, header=F, skip=1)
    names(mi) <- c(paste0("A",head))
    ind <- grep("miRNA_ID", mi[,1])
    if (length(ind)>0){
        mi <- mi[-ind,]
    }
    mi[,-1] <- sapply(mi[,-1], as.numeric)
    rownames(mi) <- mi[,1]
        
    ## CLEAN UP DATA  (REMOVE NORMAL SAMPLES; AVERAGE OR REMOVE DUPES) ## 
    mi <- straighten.up(mi, mi.map, data.type="mi")

    return(mi)
}

get.m.data <- function(countFile, fpkmFile, map, gene.info){

    head <- read.table(file = countFile, as.is=T, nrow=1)
    m <- read.table(file = countFile, as.is=T, header=F, skip=1)
    names(m) <- c("gene.id",paste0("A",head))
    rownames(m) <- gsub("\\.[0-9]+", "", m[,1])
    m$gene.id <- rownames(m)

    ## GET FPKM DATA: WILL BE USED FOR CHOOSING GENES TO KEEP BASED ON ##
    ## TPM OF 0.10 IN > 20% OF SAMPLES                                 ##
    head <- read.table(file = fpkmFile, as.is=T, nrow=1)
    f <- read.table(file = fpkmFile, as.is=T, header=F, skip=1)
    names(f) <- paste0("A",head)
    names(f)[1] <- "gene.id" 
    f$gene.id <- gsub("\\.[0-9]+", "", f$gene.id)

    ensmbl <- data.frame(ensmbl = as.character(m$gene.id), stringsAsFactors=FALSE)
    
    gene.info.in.m <- gene.info[gene.info$ensembl_gene_id %in% ensmbl[,1],]
    ## ENSEMBL ID MAPS TO MULTIPLE GENES. CONCATINATE GENES WITH AND REDUCE MULTIPLE
    ## ENTRIES OF SINGLE ENSEMBL ID TO A SINGLE ENTRY
    gene.info.in.m <- concat.dupes(gene.info.in.m)
    
    ## MERGE ENSEMBL NAMES FROM MRNA DATA WITH GENE INFORMATION FROM BIOMART
    ensmbl <- merge(ensmbl, gene.info.in.m[,c("ensembl_gene_id", "hgnc_symbol")], 
                    by.x= "ensmbl", by.y="ensembl_gene_id",
                    all.x=T, all.y=F)
    
    ## NOTE THAT ABOUT 3400 GENES DONT MAP TO A GENE NAME (NA or "")
    ## !!!!! EXCLUDING THESE GENES !!!! ##
    ## ANOTHER 27K+ ENSEMBL IDS HAVE "" FOR HGNC_SYMBOL.
    ## AT LEAST SOME ARE ANTISENSE RNA OR PROCESSED TRANSCRIPT ##
    ## ASSUMING THEY DON'T MAP FOR A REASON ##
    ind.exclude <- which(is.na(ensmbl$hgnc_symbol) | ensmbl$hgnc_symbol == "")
    
    ensmbl$hgnc_symbol[ind.exclude] <- ensmbl$ensmbl[ind.exclude]
    m <- merge(ensmbl, m, by.x = "ensmbl", by.y="gene.id", all.x = FALSE, all.y=TRUE,
               sort = FALSE)

    ## REMOVE GENES WITH NO NAMES ##
    m <- m[-grep("ENSG",m$hgnc_symbol), ]

    m <- straighten.up(m, f, map, data.type="m")
    
    return(m)
}

straighten.up <- function(data, f,  map, data.type){

    ## TRIM EXTRA COLUMNS IF DATATYPE IS MIRNA ##
    map$file_id_use = ""
    
    if (data.type == "mi"){
        names(data)[1] = "ID"
        map$file_id_use <- map$file_id_mi
    }else{
        map$file_id_use <- map$file_name_m
    }
    
    ## REMOVE NORMAL ##
    norm.ind <- grep("Normal", map$sample_type)
    norm.id <- map$file_id_use[norm.ind]

    if (length(norm.ind)>0){
        
        norm.ind <- which(gsub("^A","",names(data)) %in% norm.id)
        data <- data[, -norm.ind]
    }

    dupes <- unique(map$sample_id[duplicated(map$sample_id)])

    ## TAKE CARE OF DUPLICATE SAMPLES DEPENDING ON THE AMOUNT OF MISSINGNESS ##
    for (samp.id in dupes){

        file.ids <- map$file_id_use[map$sample_id %in% samp.id]       
        ind <- which(gsub("^A","", names(data)) %in% file.ids)

        if (length(ind) == 0){
            print(paste0("No mir data for samp.id ", samp.id, " found"))
            next
        }

        prop.zeros <- colSums(data[,ind]>0, na.rm=TRUE)/(nrow(data))

        min.zeros = min(prop.zeros)
        dif.zeros = prop.zeros - min.zeros
        rmv.dif.zeros <- which(dif.zeros > 0.1)
        
        ## REMOVE SAMPLES IF DIFFERENCE IN THE NUMBER OF ZEROS > 10% ##
        ## REMOVES SAMPLES WITH THE MOST ZEROS ##
        if (length(rmv.dif.zeros) > 0){
            data <- data[,-ind[rmv.dif.zeros]]
            ind <- ind[-rmv.dif.zeros]
        }

        ## AVERAGE EXPRESSION IF NUMBER OF REPS TWO OR MORE ##
        if (length(ind) > 1){
            data[, ind[1]] = rowMeans(data[,ind], na.rm=T)
            data <- data[,-ind[-1],]
        }       
    }

    ## LOOK FOR PROBLEMS WITH GENDER ASSIGNMENT IN MAP FILE ##
    sex.tbl <- table(map[,c("cases.0.case_id","cases.0.demographic.gender")])
    sex.tbl <- 1*(sex.tbl > 0)
    ind <- rownames(sex.tbl)[(rowSums(sex.tbl) > 1)]
    if (length(ind) > 0){
        print("Case ids ",rownames(sex.tbl)[ind], "have more than one sex assignment")
        print("Removing these samples . . .")
        ind <- which(names(countData) %in% ind)
        countData <- countData[,-ind]
    }

    ## TAKES CARE OF CASE of MULTIPLE SAMPLES AND ONE IS MISSING
    ## A GENDER ASSIGNMENT. IT WILL BE ASSIGNED THE SAME AS THE OTHER SAMPLES WITH
    ## GENEDER ASSIGNMENTS
    sex.update <- data.frame(case.id = rownames(sex.tbl),
                             gender = colnames(sex.tbl)[sex.tbl%*%c(1,2)])

    map <- merge(map, sex.update, by.x = "cases.0.case_id", by.y="case.id", all.x=T,
                 all.y=FALSE)

    ## UPDATE MAP FILE ##
    ids <- gsub("^A", "", names(data))
    if (data.type == "mi"){
        map$file_id_use <- map$file_id_mi
    }else{
        map$file_id_use <- map$file_name_m
    }

    ## REMOVE INFO FOR SAMPLES NOT IN EXPRESSION DATA ##
    map <- map[map$file_id_use %in% ids,]
    
    pref.order <- c("Primary Tumor", "Additional - New Primary",
                    "Recurrent Tumor", "Metastatic", "Additional Metastatic",
                    "Primary Blood Derived Cancer - Peripheral Blood")
    
    ## FIND DUPLICATE CASE IDS ##
    dupe.ids <- unique(map$cases.0.case_id[duplicated(map$cases.0.case_id)])
    
    if (length(dupe.ids) > 0){

        ## PICK PREFFED SAMPLE FOR EACH DUPLICATE CASE ##
        for (id in dupe.ids){
                
            inds <- grep(id, map$cases.0.case_id)
            file.ids <- map$file_id_use[inds]
            
            ind.use <- match(map$sample_type[inds], pref.order)
            ind.use <- which(ind.use == min(ind.use))
            if (length(ind.use) > 1){
                ind.use = max(ind.use)
            }

            file.id.rm <- map$file_id_use[inds[-ind.use]]

            data <- data[,-which(gsub("^A","", names(data)) %in% file.id.rm)]
        }
    }    

    ## CLEAN DATA TO EXCLUDE GENES WITH LOW COUNT INFORMATION ##
    data <- clean.data(data, f)
    
    ## REPLACE FILE ID (CURRENT HEADER) WITH CASE ID ##
    ord <- match(gsub("A","",names(data)), map$file_id_use)
    names(data)[is.na(ord)==F] <- paste0("A",map$cases.0.case_id[ord[is.na(ord)==F]])
    
    return(data)
}
    
clean.data <- function(data, f){

    ## GTEX CRITERIA IS EXCLUDE IF 20% OF SAMPLES HAVE 6 OR FEWER READS OR
    ## IF 20% OF SAMPLES HAVE TPM OF <= 0.10
    num.cols <- which(sapply(data, is.numeric))

    ## WHICH GENES HAVE 20% OR MORE OF SAMPLES WITH 6 OR FEWER READS
    rm.ind <- which( rowMeans(data[,num.cols] <=6, na.rm=T) >=.2)
    rm.id.1 <- data$ensmbl[rm.ind]

    ## DETERMINE SAMPLES WITH TPM <= .1 ##
    ## TMP_I = (FPKM_I)/SUM_J(FPKM_J)*10^6 FOR AN INDIVIDUAL
    num.cols <- which(sapply(data, is.numeric))
    f.tpm <- apply(f[,num.cols], 2, function(x) x/sum(x, na.rm=T))*10^6
    rm.ind.f <-  which(
        apply(f.tpm, 1, function(x) sum(x<=0.10, na.rm=T)>=(0.20*length(num.cols)))
        )

    genes.2.rmv <- union(f$gene.id[rm.ind.f], rm.id.1)       

    rm.ind <- which(data$ensmb %in% genes.2.rmv)
    
    data[rm.ind,-1] = NA

    return(data)
}

concat.dupes <- function(gi){

    dupe.ind <- which(duplicated(gi$ensembl_gene_id))
    dupe.names <- unique(gi$ensembl_gene_id[dupe.ind])

    for (name in dupe.names){

        ind <- which(gi$ensembl_gene_id == name)
        gi$hgnc_symbol[ind[1]] <- paste(gi$hgnc_symbol[ind], collapse=";")
        gi <- gi[-ind[-1],]
    }

    return(gi)
}

process.mi.gene <- function(data){
    
    data$miRNA = gsub("-",".",data$miRNA)
    data$miRNA = gsub("miR","mir", data$miRNA)
    data$miRNA = gsub("\\.5p|\\.3p","", data$miRNA)
    data$miRNA.alt1 <- paste0(data$miRNA, ".1")
    data$miRNA.alt2 <- paste0(data$miRNA, ".2")
    data$miRNA.alt3 <- paste0(data$miRNA, ".3")
    data$miRNA.alt4 <- paste0(data$miRNA, ".4")

    return(data)
}

est.sex <- function(m, map, gene.info){

    ## GET Y GENES THAT ARE NOT PSEUDOGENES
    ind <- which(gene.info$chromosome_name %in% c("Y"))
    genes <- gene.info$hgnc_symbol[ind]
    ind <- which(names(m) %in% genes)
    m.sex <- m[,c(1,ind)]

    map <- map[,c("cases.0.case_id","cases.0.demographic.gender")]
    names(map) = c("case.id","Sex")
    m.sex <- merge(map, m.sex, by = "case.id", all.y=T, all.x=F)
    m.sex$Sex = 1*(m.sex$Sex == "male")

    ## REMOVE GENES THAT ARE ALL NA ##
    ind <- which(colSums(is.na(m.sex))==nrow(m.sex))
    m.sex <- m.sex[,-ind]
    gene.ind <- which(names(m.sex) %in% genes)
    tmp <- lapply(gene.ind ,
                  function(x) as.vector(
                      summary(lm(m.sex[,x] ~ m.sex$Sex))$coef[2,4]))

    sig.genes <- sex.genes[which(as.numeric(sex.genes[,2])<.01), 1]
    sig.ind <- which(names(m.sex) %in% c("Sex",sig.genes))
    lm.fit <- glm(Sex ~ ., family="binomial", data = m.sex[,sig.ind])
    
    pca <- prcomp(m.sex[,sig.ind], scale=TRUE, center=TRUE)
    pca <- data.frame(pca$x[,1:10])
    pca <- data.frame(pca)
    pca$id <- m.sex$case.id
    pca$Sex <- m.sex$Sex
    summary(lm(PC1~Sex, data=pca))

    
}

get.test.pairs <- function(mi.names, m.names, mi.gene){

    m.ind <- which(mi.gene$Target.Gene %in% m.names)
    m.col <- which(names(mi.gene) == "Target.Gene")
    
    miRNA.colnames = c("miRNA",paste0("miRNA.alt",c(1:4)))
    mi.m.pairs <- c()

    for (i in 1:5){

        mi.col <- which(names(mi.gene) == miRNA.colnames[i])
        mi.ind <- which(mi.gene[,mi.col] %in% mi.names)

        ind <- intersect(mi.ind, m.ind)

        mi.m.pairs.tmp <- mi.gene[ind, c(mi.col, m.col)]
        names(mi.m.pairs.tmp) <- c("miRNA.name","gene.name")
        mi.m.pairs <- rbind(mi.m.pairs, mi.m.pairs.tmp)
    }

    ## REMOVE DUPES ##
    nms <- paste(mi.m.pairs$miRNA.name, mi.m.pairs$gene.name, sep=".")
    dupe.ind <- which(duplicated(nms))
    mi.m.pairs <- mi.m.pairs[-dupe.ind,]

    return(mi.m.pairs)            
}

    
test.pairs <- function(mi, m, mi.m.pairs){

    out <- c()

    ## PROCEDURE:
    ## 1)LOG TRANSFORM BOTH MRNA AND MIRNA
    ## 2)PERFORM PCA ON MRNA AND MIRNA
    ## 3)CREATE ADJUST MIRNA DATASET THAT REGRESSES OUT PCS FROM MIRNA
    ## 4)REGRESS LOG-MRNA ~ MRNA.PC1. . .+ MIRNA.ADJUSTED + SEX + MIRNA.ADJUSTED*SEX
    ## 5)REGRESS LOG-MRNA ~ MRNA.PC1. . .+ MIRNA.ADJUSTED.1 . . .+ SEX +
    ## . . . SEX*MIRNA.ADJUSTED.1+ . . .
    ## 6)USE MIRNA AS DICOTOMOUS (EXPRESSED/NOT EXPRESSED) AND REPEAT ABOVE

    ## MERGE MI AND M
    data <- merge(m, mi, by = "case.id", all = FALSE)

    ## REMOVE COLUMNS WITH ALL NAs FROM DATA##
    ind <- which(colSums(is.na(data)) == nrow(data))
    if (length(ind) != 0){
        rm.nms <- names(data)[ind]
        data = data[,-ind]
    }

    ## REMOVE REMOVED MIRNAS AND GENES FROM MI.M.PAIRS
    ind <- which(mi.m.pairs$gene.name %in% rm.nms |
                 mi.m.pairs$miRNA.name %in% rm.nms)
    if (length(ind) != 0){
        mi.m.pairs <- mi.m.pairs[-ind,]
    }
    
    ## REDUCE DATA TO ONLY HAVE NECESSARY GENES AND MIRNAS ##
    genes <- unique(mi.m.pairs$gene.name)
    mirs <- unique(mi.m.pairs$miRNA.name)
    ind <- which(names(data) %in% c(genes, mirs, "Sex"))
    data <- data[,c(1,ind)]

    m.names.ind <- which(names(data) %in% genes)
    m.names <- names(data)[ m.names.ind ]
    mi.names.ind <- which(names(data) %in% mirs)
    mi.names <- names(data)[ mi.names.ind ]
    
    ## LOG TRANSFOR MIRNA AND EXPRESSION LEVELS ##
    sex.ind <- which(names(data) == "Sex")
    data.log <- data
    tmp <- data.log[,-c(1,sex.ind)] + 1
    tmp <- log(tmp)
    data.log[,-c(1,sex.ind)] <- tmp       
    
    ## PCA ##
    pca.m <- prcomp(data.log[,m.names.ind], scale=TRUE, center=TRUE)
    pca.m <- data.frame(pca.m$x[,1:10])
    pca.m$case.id <- data.log$case.id

    ## ADD MRNA PCS TO DATA ##    
    data <- merge(data, pca.m, by="case.id")
    
    ## ADD MIRESIDS TO DATA ##

    data$Sex <- 1*(data$Sex == "male")
    
    out <- c()
    #out.mires <- c()
    out.gene <- c()
    gene.count <- 1
    for (gene in genes){

        if (gene.count%%25 == 0){ print(paste("Working on gene", gene.count)) }
        gene.count = gene.count+1
        gene.ind <- which(names(data) == gene)
        ## CREATE DATASET WITH ONLY ONE GENE AND ASSOCIATED MIRNAS
        mir.names <- mi.m.pairs$miRNA.name[mi.m.pairs$gene.name == gene]
        mir.ind <- which(names(data) %in% mir.names)
        ## mir.res.ind <- which(names(data) %in% paste(mir.names, "res", sep="."))
        
        lms <- lapply(mir.ind,
                      function(x) as.vector(t(
                          summary(lm(data[,gene.ind] ~ data$PC1 +
                                     data$PC2 + data$PC3 + data$PC4 + data$PC5 +
                                     data$Sex + data[,x]))$coef[7:8,])))
        lms.mtx <- do.call(rbind, lms)
        
        lms.int <- lapply(mir.ind,
                          function(x) as.vector(t(
                              summary(lm(data[,gene.ind] ~ data$PC1 +
                                      data$PC2 + data$PC3 + data$PC4 + data$PC5 +
                                      data$Sex + data[,x] +
                                         data$Sex*data[,x]))$coef[c(7,9),])))

        lms.int.mtx <- do.call(rbind, lms.int)

        out <- rbind(out, cbind(gene, names(data)[mir.ind], lms.mtx, lms.int.mtx))

        ## TEST ALL MIRS FOR A GENE AT ONCE ##
        if (length(mir.ind) > 1){
            
            models <- build.model(mir.ind, gene.ind)
            lm.red <- eval(parse(text=models[[1]]))
            lm.full <- eval(parse(text=models[[2]]))
            anva <- anova(lm.red, lm.full)
            
            out.gene <- rbind(out.gene, c(gene, length(mir.ind), anva$P[2]))
        }else{
            out.gene <- rbind(out.gene, c(gene, 1, lms.int[8]))
        }
        
        ## REPEAT WITH THE RESIDUALIZED MIR DATA ##
        ## lms <- lapply(mir.res.ind,
        ##             function(x) as.vector(t(
        ##                summary(lm(data[,gene.ind] ~ data$PC1 +
        ##                            data$PC2 + data$PC3 + data$PC4 + data$PC5 +
        ##                           data$Sex + data[,x]))$coef[7:8,])))
        ## lms.mtx <- do.call(rbind, lms)
        ## lms.int <- lapply(mir.res.ind,
        ##                  function(x) as.vector(t(
        ##                      summary(lm(data[,gene.ind] ~ data$PC1 +
        ##                              data$PC2 + data$PC3 + data$PC4 + data$PC5 +
        ##                              data$Sex + data[,x] +
        ##                              data$Sex*data[,x]))$coef[c(7,9),])))
        ## lms.int.mtx <- do.call(rbind, lms.int)
        ## out.mires <- rbind(out.mires, cbind(gene, names(data)[mir.res.ind], lms.mtx, lms.int.mtx))

      
    }
    
    out <- data.frame(out)
    names(out)[1:2] = c("gene","mir")
    names(out)[-c(1:2)] <- as.vector(do.call(cbind,
        lapply(list("miRNA","Sex","Sex.g.Int","Int"),
               function(x) paste(x,c("Est","SE","t","P"),
                                 sep="."))))

    out$Int.fdr <- p.adjust(as.numeric(as.character(out$Int.P)), method="BY")

    ## out.mires <- data.frame(out.mires)
    ## names(out.mires)[1:2] = c("gene","mir")
    ## names(out.mires)[-c(1:2)] <- as.vector(do.call(cbind,
    ##     lapply(list("miRNA","Sex","Sex.g.Int","Int"),
    ##           function(x) paste(x,c("Est","SE","t","P"),
    ##                             sep="."))))

    ## out.mires$Int.fdr <-  p.adjust(as.numeric(as.character(out.mires$Int.P)),
    ##                              method="BY")
    

    return(list(out, out.gene))
}

build.model <- function(mir.res.ind, gene.ind){

    model <- paste("lm(data[,",gene.ind,"]~ data$PC1 + data$PC2 +",
                   "data$PC3 + data$PC4 + data$PC5 +data$Sex +")
    datas <- sapply(mir.res.ind, function(x) paste("data[,",x,"]", sep=""))
    model <- paste(model, paste(datas, collapse=" + "))
    model.red <- paste(model, ")")
    ints <-  paste("data$Sex", datas, sep="*")
    model.full <- paste(model, " + ", paste(ints, collapse=" + "), ")")

    return(list(model.red = model.red, model.full = model.full))
    
}

## test.pairs(mi.cl, m.cl, mi.m.pairs)

di.test.pairs <- function(mi, m, mi.m.pairs){

    out <- c()

    ## PROCEDURE:
    ## 1)LOG TRANSFORM BOTH MRNA AND MIRNA
    ## 2)PERFORM PCA ON MRNA AND MIRNA
    ## 3)CREATE ADJUST MIRNA DATASET THAT REGRESSES OUT PCS FROM MIRNA
    ## 4)REGRESS LOG-MRNA ~ MRNA.PC1. . .+ MIRNA.ADJUSTED + SEX + MIRNA.ADJUSTED*SEX
    ## 5)REGRESS LOG-MRNA ~ MRNA.PC1. . .+ MIRNA.ADJUSTED.1 . . .+ SEX +
    ## . . . SEX*MIRNA.ADJUSTED.1+ . . .
    ## 6)USE MIRNA AS DICOTOMOUS (EXPRESSED/NOT EXPRESSED) AND REPEAT ABOVE

    ## MERGE MI AND M
    data <- merge(m, mi, by = "case.id", all = FALSE)

    ## REMOVE COLUMNS WITH ALL NAs FROM DATA##
    ind <- which(colSums(is.na(data)) == nrow(data))
    if (length(ind) != 0){
        rm.nms <- names(data)[ind]
        data = data[,-ind]
    }

    ## REMOVE REMOVED MIRNAS AND GENES FROM MI.M.PAIRS
    ind <- which(mi.m.pairs$gene.name %in% rm.nms |
                 mi.m.pairs$miRNA.name %in% rm.nms)
    if (length(ind) != 0){
        mi.m.pairs <- mi.m.pairs[-ind,]
    }
    
    ## REDUCE DATA TO ONLY HAVE NECESSARY GENES AND MIRNAS ##
    genes <- unique(mi.m.pairs$gene.name)
    mirs <- unique(mi.m.pairs$miRNA.name)
    ind <- which(names(data) %in% c(genes, mirs, "Sex"))
    data <- data[,c(1,ind)]

    m.names.ind <- which(names(data) %in% genes)
    m.names <- names(data)[ m.names.ind ]
    mi.names.ind <- which(names(data) %in% mirs)
    mi.names <- names(data)[ mi.names.ind ]
    
    ## LOG TRANSFOR MIRNA AND EXPRESSION LEVELS ##
    sex.ind <- which(names(data) == "Sex")
    data.log <- data
    tmp <- data.log[,-c(1,sex.ind)] + 1
    tmp <- log(tmp)
    data.log[,-c(1,sex.ind)] <- tmp       
    
    ## PCA ##
    pca.m <- prcomp(data.log[,m.names.ind], scale=TRUE, center=TRUE)
    pca.m <- data.frame(pca.m$x[,1:10])
    pca.m$case.id <- data.log$case.id

    ## ADD MRNA PCS TO DATA ##    
    data <- merge(data, pca.m, by="case.id")
    
    ## CODE SEX AS 0=FEMALE, 1= MALE
    data$Sex <- 1*(data$Sex == "male")

    ## TURN MIRNA TO DICOTOMOUS DATA ##
    noexpress.cutoff = 10
    ind <- names(data) %in% mirs 
    data[,ind] <- 1*(data[,ind] > noexpress.cutoff)
    ## REMOVE MIRS THAT DONT HAVE AT LEAST 10 EXPRESSED OR UNEXPRESSED SAMPLES ##
    n.expressed <- colSums(data[,ind])
    rmv.col <- which(n.expressed < 10 | n.expressed > (nrow(data)-10))
    data <- data[, -which(ind)[rmv.col]]

    
    out <- c()
    ## out.mires <- c()
    out.gene <- c()
    gene.count <- 1
    for (gene in genes){

        if (gene.count%%25 == 0){ print(paste("Working on gene", gene.count)) }
        gene.count = gene.count+1
        gene.ind <- which(names(data) == gene)
        ## CREATE DATASET WITH ONLY ONE GENE AND ASSOCIATED MIRNAS
        mir.names <- mi.m.pairs$miRNA.name[mi.m.pairs$gene.name == gene]
        mir.ind <- which(names(data) %in% mir.names)

        ## IT MAY BE THAT NONE OF A GENES MIRS HAD ENOUGH SUBJECTS EXPRESSED
        ## AND UNEXPRESSED SUBJECTS TO TEST GENES
        if (length(mir.ind) == 0){ next }
        ## mir.res.ind <- which(names(data) %in% paste(mir.names, "res", sep="."))
        
        lms <- lapply(mir.ind,
                      function(x) as.vector(t(
                          summary(lm(data[,gene.ind] ~ data$PC1 +
                                     data$PC2 + data$PC3 + data$PC4 + data$PC5 +
                                     data$Sex + data[,x]))$coef[7:8,])))
        lms.mtx <- do.call(rbind, lms)
        
        lms.int <- lapply(mir.ind,
                          function(x) as.vector(t(
                              summary(lm(data[,gene.ind] ~ data$PC1 +
                                      data$PC2 + data$PC3 + data$PC4 + data$PC5 +
                                      data$Sex + data[,x] +
                                         data$Sex*data[,x]))$coef[c(7,9),])))

        lms.int.mtx <- do.call(rbind, lms.int)

        out <- rbind(out, cbind(gene, names(data)[mir.ind], lms.mtx, lms.int.mtx))

        ## TEST ALL MIRS FOR A GENE AT ONCE ##
        if (length(mir.ind) > 1){
            
            models <- build.model(mir.ind, gene.ind)
            lm.red <- eval(parse(text=models[[1]]))
            lm.full <- eval(parse(text=models[[2]]))
            anva <- anova(lm.red, lm.full)
            
            out.gene <- rbind(out.gene, c(gene, length(mir.ind), anva$P[2]))
        }else{
            out.gene <- rbind(out.gene, c(gene, 1, lms.int[8]))
        }                    
    }
    
    out <- data.frame(out)
    names(out)[1:2] = c("gene","mir")
    names(out)[-c(1:2)] <- as.vector(do.call(cbind,
        lapply(list("miRNA","Sex","Sex.g.Int","Int"),
               function(x) paste(x,c("Est","SE","t","P"),
                                 sep="."))))

    out$Int.fdr <- p.adjust(as.numeric(as.character(out$Int.P)), method="BY")

    return(list(out, out.gene))
}
