## ---------------------------
##
## Script name: Preprocess blast and bakta results file
##
## Purpose of script:
##
##
## Author: Sneha Sundar
##
## 
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


#################################################
# PACKAGES AND SOURCE FUNCTIONS

library(tidyverse)
library(here)
source(here("code/000.tools/helpers.R"))


##################################################
# DEFINE GLOBAL VARIABLES

#Ordered F plasmid tra region genes (koraiman2018)
TRA_GENES <- c("tram","traj","tray","traa","tral","trae","trak","trab","trap","trbd","trbg","trav","trar","trac","trbi","traw","trau","trbc","tran","trbe","traf","trba","arta","traq","trbb","trbj","trbf","trah","trag","tras","trat","trad","trai","trax","fino")

#tra gene metadata file (query metadata)
TRA_GENE_METADATA_FP <- here("data/processed/002.find_genes/reference_tra_genes/tra_ref_genes_metadata_nodup_manual_35_genes.tsv")

DATASET <-  "pires1200/ecoli"

BLAST_RESULTS_FOLDER <-  here(paste0("data/processed/002.find_genes/",DATASET,"/blast_results"))
BAKTA_RESULTS_FOLDER <- here(paste0("data/processed/002.find_genes/",DATASET,"/bakta_results"))

#read in list of genomes to analyse
GENOMES <- read_lines(here("data/raw/pires_1200/BioSample_ecoli.txt"))

#how many best hits to keep and at what level
N_TOP_HITS <- 1
LEVEL <- "contig"

#remove hits with pid < PID_THRES
PID_THRES <- 50

#remove hits with cov < PID_THRES
COV_THRES <- 30

#outfolder
OUTFOLDER <-  here(paste0("data/processed/002.find_genes/",DATASET))

###########################################


for(genome in GENOMES[1:1200]){
  
  ################ PREPROCESS BLAST RESULTS ###################
  
  blast_res <- load_raw_blast6out(BLAST_RESULTS_FOLDER, genome, extension = ".tblastn6out")
  
  tra_gene_metadata <- load_tra_gene_metadata(TRA_GENE_METADATA_FP)
  
  #add tra gene metadata
  blast_res <- add_tra_gene_metadata(blast_res,tra_gene_metadata)
  
  #compute query coverage
  blast_res <- compute_qcov(blast_res)
  
  #pick only the best n_hits 
  blast_res <- pick_top_n_hits(blast_res, on = LEVEL, n = N_TOP_HITS)
  
  #keep only actual tra genes (removing vagd, orf1 and proq/fino)
  blast_res <- blast_res %>% filter(gene_name %in% TRA_GENES)
  
  #plot pid and qcov of hits as a sanity check
  p <- plot_pid_vs_qcov(blast_res, PID_THRES, COV_THRES)
  pdf(paste0(OUTFOLDER,"/pid_vs_qcov_unfiltered/",genome,".pdf"))
  print(p)
  dev.off()
  
  #remove low pid and qcov hits
  blast_res <- keep_good_hits(blast_res, PID_THRES, COV_THRES)
  
  if(nrow(blast_res) != 0){
  #assign strand
  strand_assign_df <- t(mapply(assign_strand,blast_res$tstart,blast_res$tend))
  blast_res$tstart <- as.integer(strand_assign_df[,1])
  blast_res$tend <- as.integer(strand_assign_df[,2])
  blast_res["strand"] <- strand_assign_df[,3]
  #calculate target seq protein length
  blast_res["protein_len"] <- (abs(blast_res$tend-blast_res$tstart)+1)/3
  write.table(blast_res,file = paste0(OUTFOLDER,"/preprocessed_blast_res/bitscore_qcov_pid/",genome,".tblastn6outmod"),quote=F,sep='\t',row.names=F,col.names = F)
  }else{
    write_lines(genome,paste0(OUTFOLDER,"/preprocessed_blast_res/bitscore_qcov_pid/","empty.txt"),append = TRUE)
  }
  
  
  
  ################ PREPROCESS BAKTA RESULTS TO INCLUDE ONLY CONTIGS THAT HAVE F_LIKE TRA GENES IDENTIFIED BY BLAST ###################
  
  if(nrow(blast_res)!=0){
  bakta_res <- load_raw_bakta(BAKTA_RESULTS_FOLDER,genome, extension = ".tsv")
  
  bakta_res <- bakta_res %>%  mutate(gene = str_to_lower(na_if(gene,"")))
  
  contigs_with_tra_hits <- get_tra_gene_contigs(blast_res)
  
  
  bakta_res <- get_bakta_contigs(bakta_res,contigs_with_tra_hits)
  
  #add genome info 
  bakta_res <- cbind(genome, bakta_res)
  
  if(nrow(bakta_res!=0)){
  write.table(bakta_res,file = paste0(OUTFOLDER,"/preprocessed_bakta_res/bitscore_qcov_pid/",genome,".baktamod"),quote=F,sep='\t',row.names=F,col.names = F)
  }else{
    #only those genomes with tra hits but no bakta hits (weird case)
    write_lines(genome,paste0(OUTFOLDER,"/preprocessed_bakta_res/bitscore_qcov_pid/","empty.txt"),append = TRUE)
  }
  
}
}
write_lines(colnames(blast_res),paste0(OUTFOLDER,"/preprocessed_blast_res/","colnames.txt"))
write_lines(colnames(bakta_res),paste0(OUTFOLDER,"/preprocessed_bakta_res/","colnames.txt"))
