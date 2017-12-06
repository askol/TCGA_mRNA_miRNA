setwd("/group/stranger-lab/askol/TCGA/Expression_Analysis/")

projects = read.table(file = "../TCGA_Target_projects.txt", header=TRUE, as.is=TRUE, sep="\t")
projects = projects[grep("TCGA",projects[,1]),1]

job.files = c()
for (project in projects){

    job.file = paste0("DE_analysis_",project,".pbs")
    job.files = c(job.files, job.file)

    write(file = job.file, "#!/bin/bash")
    write(file = job.file, paste("#PBS -l nodes=1:ppn=1,mem=16gb"), append=TRUE)
    write(file=job.file, "#PBS -l walltime=36:00:00", append=TRUE)
    write(file = job.file, "#PBS -o /group/stranger-lab/askol/TCGA/Expression_Analysis", append=TRUE)
    write(file = job.file, "#PBS -j oe", append=TRUE)
    write(file = job.file, paste0("#PBS -N DE_Analysis_",project), append=TRUE)
    write(file = job.file, paste0("R CMD BATCH  '--args ",project,
              "' /group/stranger-lab/askol/TCGA/Code/process_expression_data.r ",
              "/group/stranger-lab/askol/TCGA/Expression_Analysis/",project,".out"),
              append=TRUE)
    
    system(paste("chmod +x", job.file))
}

for (job in job.files){
    system(paste("qsub",job))
}
 
