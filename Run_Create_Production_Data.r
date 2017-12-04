projects = read.table(file = "/group/stranger-lab/askol/TCGA/TCGA_Target_projects.txt", header=TRUE, as.is=TRUE, sep="\t")
projects = projects[grep("TCGA",projects[,1]),1]

job.files = c()
setwd("/group/stranger-lab/askol/TCGA/Make_Production_Data")

for (project in projects){

    skip.cancers <- c("TCGA-BRCA","TCGA-OV", "TCGA-CESC", "TCGA-PRAD", "TCGA-TGCT", "TCGA-UCS", "TCGA-UCEC")
    if (project %in% skip.cancers){ next }
    job.file = paste0("create_production_data",project,".pbs")
    job.files = c(job.files, job.file)

    write(file = job.file, "#!/bin/bash")
    write(file = job.file, paste("#PBS -l nodes=1:ppn=1,mem=8gb"), append=TRUE)
    write(file=job.file, "#PBS -l walltime=36:00:00", append=TRUE)
    write(file = job.file, "#PBS -o /group/stranger-lab/askol/TCGA/Make_Production_Data", append=TRUE)
    write(file = job.file, "#PBS -j oe", append=TRUE)
    write(file = job.file, paste0("#PBS -N Create_Prod_Data_",project), append=TRUE)
    write(file = job.file, paste0("R CMD BATCH  '--args ",project,
              "' /group/stranger-lab/askol/TCGA/Code/Create_Production_Data.r ",
              "/group/stranger-lab/askol/TCGA/Make_Production_Data/",project,".out"),
              append=TRUE)
    
    system(paste("chmod +x", job.file))
}

for (job in job.files){
    system(paste("qsub",job))
}
 
