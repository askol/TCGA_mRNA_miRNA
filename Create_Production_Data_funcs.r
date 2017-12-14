make.mi.map <- function(project){

    METADIR <- "/group/stranger-lab/askol/TCGA/TCGA_MetaData"
      
    meta.miRNA.file <- paste0(METADIR,"/",project,"_miRNA_meta.txt")
      
    mi.map <- read.table(file=meta.miRNA.file, header=T, as.is=T, sep="\t")
    names(mi.map) <- gsub("file_id", "file_id_mi", names(mi.map))
    names(mi.map) <- gsub("cases.0.samples.0.","", names(mi.map))

    if (!any(grep("gender", names(mi.map)))){
        mi.map$cases.0.demographic.gender = NA
    }
    
    mi.map <- gender.fill(mi.map)
    
    return(mi.map)
}

make.m.map <- function(project){

     METADIR <- "/group/stranger-lab/askol/TCGA/TCGA_MetaData"
     meta.mRNA.file <- paste0(METADIR,"/",project,"_mRNA_meta.txt")
     
     m.map <- read.table(file=meta.mRNA.file, header=T, as.is=T, sep="\t")
     
     names(m.map) <- gsub("file_name", "file_name_m", names(m.map))
     m.map$file_name_m = gsub("\\..*","",m.map$file_name_m)
     
     names(m.map) <- gsub("cases.0.samples.0.","", names(m.map))

     m.map <- gender.fill(m.map)
     
     return(m.map)
 }

gender.fill <- function(map){

    ## LOOK FOR PROBLEMS WITH GENDER ASSIGNMENT IN MAP FILE ##
    sex.tbl <- table(map[,c("cases.0.case_id","cases.0.demographic.gender")])
    ## IF SEX IS MISSING, THEN THERE WILL BE A COLUMN WITH NO COLNAME THAT HAS
    ## A 1 WHERE NO MALE OR FEMALE VALUE WAS STATED. REMOVE THIS COLUMN SINCE THE
    ## MALE AND FEMALE COLUMNS WILL EACH HAVE 0
    ind <- which(colnames(sex.tbl) == "")
    if (length(ind) >0){
        sex.tbl = sex.tbl[, -ind]
    }
    sex.tbl <- 1*(sex.tbl > 0)
    ind <- rownames(sex.tbl)[(rowSums(sex.tbl) > 1)]
    if (length(ind) > 0){
        print("Case ids ",rownames(sex.tbl)[ind], "have more than one sex assignment")
        print("Removing these samples . . .")       
    }

    ## TAKES CARE OF CASE of MULTIPLE SAMPLES AND ONE IS MISSING
    ## A GENDER ASSIGNMENT. IT WILL BE ASSIGNED THE SAME AS THE OTHER SAMPLES WITH
    ## GENEDER ASSIGNMENTS
    sex.update <- data.frame(case.id = rownames(sex.tbl),
                             gender = NA)
    ind = sex.tbl%*%c(1,2)
    sex.update$gender[ind !=0] <- colnames(sex.tbl)[ind[ind!=0]]

    map <- merge(map, sex.update, by.x = "cases.0.case_id", by.y="case.id", all.x=T,
                 all.y=FALSE)

    return(map)
}

make.map.files <- function(project){
    
    m.map <- make.mi.map.file(project)
    mi.map <- make.m.map.file(project)
    
    return(list(m.map = m.map, mi.map = mi.map))

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
    mi <- straighten.up(mi, "", mi.map, data.type="mi")

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
    head <- read.table(file = fpkmFile, as.is=T, nrow=2)
    skip = 1
    ## A COUPLE FILES HAVE TWO HEADER LINES, REMOVE IF NECESSARY ##
    if (grep("ENSG", head[2,1]) == T){
        skip=2
    }
    head <- head[1,]
    
    f <- read.table(file = fpkmFile, as.is=T, header=F, skip=skip)

    names(f) <- paste0("A",head)
    names(f)[1] <- "gene.id" 
    f$gene.id <- gsub("\\.[0-9]+", "", f$gene.id)

    ## ensmbl <- data.frame(ensmbl = as.character(m$gene.id), stringsAsFactors=FALSE)
    
    ## gene.info.in.m <- gene.info[gene.info$ensembl_gene_id %in% ensmbl[,1],]

    ## ENSEMBL ID MAPS TO MULTIPLE GENES. CONCATINATE GENES WITH AND REDUCE MULTIPLE
    ## ENTRIES OF SINGLE ENSEMBL ID TO A SINGLE ENTRY
    ## gene.info.in.m <- concat.dupes(gene.info.in.m)
    
    ## MERGE ENSEMBL NAMES FROM MRNA DATA WITH GENE INFORMATION FROM BIOMART
    ## ensmbl <- merge(ensmbl, gene.info.in.m[,c("ensembl_gene_id", "hgnc_symbol")], 
    ##                by.x= "ensmbl", by.y="ensembl_gene_id",
    ##                all.x=T, all.y=F)
    
    ## NOTE THAT ABOUT 3400 GENES DONT MAP TO A GENE NAME (NA or "")
    ## !!!!! EXCLUDING THESE GENES !!!! ##
    ## ANOTHER 27K+ ENSEMBL IDS HAVE "" FOR HGNC_SYMBOL.
    ## AT LEAST SOME ARE ANTISENSE RNA OR PROCESSED TRANSCRIPT ##
    ## ASSUMING THEY DON'T MAP FOR A REASON ##

    ## ind.exclude <- which(is.na(ensmbl$hgnc_symbol) | ensmbl$hgnc_symbol == "")
    
    ## ensmbl$hgnc_symbol[ind.exclude] <- ensmbl$ensmbl[ind.exclude]
    ## m <- merge(ensmbl, m, by.x = "ensmbl", by.y="gene.id", all.x = FALSE, all.y=TRUE,
    ##           sort = FALSE)

    ## REMOVE GENES WITH NO NAMES ##
    ## m <- m[-grep("ENSG",m$hgnc_symbol), ]

    m <- straighten.up(m, f, map, data.type="m")
    
    return(m)
}

straighten.up <- function(data, f="",  map, data.type){

    ## TRIM EXTRA COLUMNS IF DATATYPE IS MIRNA ##
    map$file_id_use = ""
    
    if (data.type == "mi"){
        names(data)[1] = "ID"
        map$file_id_use <- map$file_id_mi
    }else{
        names(data)[1] <- "ensmbl"
        map$file_id_use <- map$file_name_m
        
    }
    
    ## REMOVE NORMAL ##
    norm.ind <- grep("Normal", map$sample_type)
    norm.id <- map$file_id_use[norm.ind]

    if (length(norm.ind)>0){
        
        norm.ind <- which(gsub("^A","",names(data)) %in% norm.id)
        data <- data[, -norm.ind]
        print(paste0("Removing ",length(norm.ind), " normal samples."))
    }

    dupe.inds <- which(duplicated(map$sample_id))
    dupes <- unique(map$sample_id[dupe.inds])
    if (length(dupes)>0){
        print(paste0(length(dupe.inds)," duplicates from ",length(dupes), " samples"))
    }
    
    ## TAKE CARE OF DUPLICATE SAMPLES DEPENDING ON THE AMOUNT OF MISSINGNESS ##
    for (samp.id in dupes){

        file.ids <- map$file_id_use[map$sample_id %in% samp.id]       
        ind <- which(gsub("^A","", names(data)) %in% file.ids)

        if (length(ind) == 0){
            print(paste0("No mir/mrna data for samp.id ", samp.id, " found"))
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
    dupe.ind <- which(duplicated(map$cases.0.case_id))
    dupe.ids <- unique(map$cases.0.case_id[dupe.ind])

    if (length(dupe.ids) > 0){

        print(paste0(length(dupe.ind), " duplicates by case, made up of ",
                     length(dupe.ids), " cases."))

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
    if (data.type == "mi"){
        data <- clean.data(data, f="", data.type=data.type)
    }else{
        data <- clean.data(data, f=f, data.type=data.type)
    }

    ## REPLACE FILE ID (CURRENT HEADER) WITH CASE ID ##
    ord <- match(gsub("A","",names(data)), map$file_id_use)
    names(data)[is.na(ord)==F] <- paste0("A",map$cases.0.case_id[ord[is.na(ord)==F]])

    ## CONSOLIDATE MAP TO HAVE ONLY A SINGLE ENTRY PER CASE.ID ##
    map <- map[map$cases.0.case_id %in% gsub("A","",names(data)),]
    tbl <-  table(map$cases.0.demographic.gender)
    print(paste("Number", paste0(names(tbl),"s : ", tbl), collapse="; "))

    ## REMOVE CASES WITH NO GENDER ##
    ind <- map$cases.0.case_id[is.na(map$cases.0.demographic.gender)]
    ind <- which(names(data) %in% ind)
    if (length(ind) > 0){
        data <- data[,-ind]
    }
        
    return(data)
}
    
clean.data <- function(data, f="", data.type = "m"){

  
    num.cols <- which(sapply(data, is.numeric))
    rna.id.col <- which(names(data) %in% c("ID","ensmbl"))
    
    ## MIRNA DOES NOT USE TPM CUTOFF BECAUSE THE LENGTH OF MIRNA ARE LESS
    ## THAN THE READ LENGTH (DON'T EXPECT NUMBER OF READS TO INCREASE WITH MIRNA
    ## LENGTH
    if (data.type == "mi"){
        rm.ind <- which( rowMeans(data[,num.cols] <=5, na.rm=T) >=.8)
        rm.ind.1 <- which( apply(data[, num.cols], 1, sd, na.rm=T) < 4)
        rm.id <- data[intersect(rm.ind, rm.ind.1), rna.id.col]
        rm.ind <- which(data[, rna.id.col] %in% rm.id)
    }

    ## GTEX CRITERIA IS EXCLUDE IF 20% OF SAMPLES HAVE 6 OR FEWER READS OR
    ## IF 20% OF SAMPLES HAVE TPM OF <= 0.10 {CHANGED TO 80% BECAUSE GOAL IS
    ## IS DIFFERENTIAL EXPRESSION}
    if (data.type == "m"){
    
        ## WHICH GENES HAVE 20% OR MORE OF SAMPLES WITH 6 OR FEWER READS
        rm.ind <- which( rowMeans(data[,num.cols] <=6, na.rm=T) >=.8)
        rm.id.1 <- data[rm.ind, rna.id.col]

        ## DETERMINE SAMPLES WITH TPM <= .1 ##
        ## TMP_I = (FPKM_I)/SUM_J(FPKM_J)*10^6 FOR AN INDIVIDUAL
        num.cols <- which(sapply(f, is.numeric))
        f.tpm <- apply(f[,num.cols], 2, function(x) x/sum(x, na.rm=T))*10^6
        rm.ind.f <-  which(
            apply(f.tpm, 1, function(x) sum(x<=0.10, na.rm=T)>=(0.80*length(num.cols)))
            )

        rm.id <- union(f$gene.id[rm.ind.f], rm.id.1)
        rm.ind <- which(data[,rna.id.col] %in% rm.id)
        
    }
    
    print(paste0("Setting ",length(rm.ind), " of ",nrow(data), " to missing. ",
                 nrow(data)-length(rm.ind), " are not missing"))
    
    data <- data[-rm.ind , ]
    
    ## EXCLUDE GENE/MIRNA IF 10 or FEWER ENTRIES HAVE DATA (E.G. REST NAS)
    rm.genes <- which(rowSums(!is.na(data[,num.cols]))<=10)
    if (length(rm.genes) >0){
        data <- data[-rm.genes,]
        print(paste0("Removing ",length(rm.genes), " genes due to too many NAs. ",
                     "Fewer than 10 subjects with counts. "))
        print(paste0(length(num.cols) - length(rm.genes), " genes/mirs remain"))
    }
    
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

get.xy.mirs <- function(){

    ## GET MIR LOCATIONS ##
    mir.locs <- read.table(file = "/group/stranger-lab/askol/TCGA/hsa.gff3", skip = 13,
                           header=FALSE, sep="\t")
    mir.locs[,9] = gsub(".*Name=","",mir.locs[,9])
    mir.locs[,9] = gsub(";.*","", mir.locs[,9])
    mir.locs <- mir.locs[,c(9,1)]
    names(mir.locs) <- c("mir","chr")
    xy.mirs <- mir.locs$mir[grep("x|y|X|Y",mir.locs$chr)]
    xy.mirs <- unique(xy.mirs)
    xy.mirs <- gsub("-","\\.",xy.mirs)

    ## REMOVE 3P/5P FROM END OF MIRNA NAME ##
    xy.mirs <- gsub("miR","mir", xy.mirs)
    xy.mirs <- gsub("\\.3p$|\\.5p$","", xy.mirs)

    return(xy.mirs)
}

get.mir.info <- function(){

    ## GET MIR LOCATIONS FOR PRIMARY TRANSCRIPTS ##
    mir.locs <- read.table(file = "/group/stranger-lab/askol/TCGA/hsa.gff3", skip = 13,
                           header=FALSE, sep="\t")
    ## KEEP ONLY PRIMARY TRANSCRIPT ENTRIES ##
    keep.ind <- grep("primary",mir.locs[,3])
    mir.locs <- mir.locs[ keep.ind, c(9,1,4,5)]
    names(mir.locs) <- c("mir","chr","start","end")
    mir.locs$mir = gsub(".*Name=","",mir.locs$mir)
    mir.locs$mir = gsub(";.*","", mir.locs$mir)
    
    return(mir.locs)
}

get.mir.info.trim.name <- function(){

    mir.locs <- get.mir.info()

    ##    
    mirs <- strsplit(split="-", mir.locs$mir)
    mirs <- sapply(mirs, function(x){paste(x[1:3], collapse="-")})
    mir.locs$mir.long <- mir.locs$mir
    mir.locs$mir <- mirs

    chrs <- with(mir.locs, tapply(chr, mir, paste, collapse=","))
    chrs <- data.frame(mir = rownames(chrs), chrs=chrs)
    
    starts <- with(mir.locs, tapply(start, mir, paste, collapse=","))
    starts <- data.frame(mir = rownames(starts), starts=starts)

    ends <- with(mir.locs, tapply(end, mir, paste, collapse=","))
    ends <- data.frame(mir = rownames(ends), ends=ends)

    mir.longs <- with(mir.locs, tapply(mir.long, mir, paste, collapse=","))
    mir.longs <- data.frame(mir = rownames(mir.longs), ends=ends)
    
    mir.locs <- merge(mir.locs, chrs, by="mir", all.x=T, all.y=F)
    mir.locs <- merge(mir.locs, starts, by="mir", all.x=T, all.y=F)
    mir.locs <- merge(mir.locs, ends, by="mir", all.x=T, all.y=F)
    mir.locs <- merge(mir.locs, mir.longs, by="mir", all.x=T, all.y=F)
    
    ## REMOVE DUPES ##
    dupe.ind <- which(duplicated(mir.locs$mir))
    mir.locs <- mir.locs[-dupe.ind , -which(names(mir.locs) %in%
                                            c("chr","start","end"))]
    
    return(mir.locs)
}
    
