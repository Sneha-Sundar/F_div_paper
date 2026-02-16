## ---------------------------
##
## Script name: Create blast summary tables
##
## Purpose of script:
##
## Filter preprocessed blast tables based on PID and QCOV thresholds
## Remove lone TrbH and traT containting contigs
## Make contig-level and genome level summary tables describing number of tra gene hits, median pid, median qcov and metadata on the genomes.
##
## Author: Sneha Sundar
##
## Date Created: 2024-02-02
##
## ---------------------------
##
## Notes: The input preprocessed tables contain top hits with min 50% pid and min 30% qcov.
##
##
## ---------------------------


#################################################
# PACKAGES AND SOURCE FUNCTIONS


library(here)
library(tidyverse)
source(here("code/000.tools/plot.R"))


#################################################
# INPUT


INPUT_BLAST_FILE <-
  here(
    "data/processed/002.find_genes/pires1200/ecoli/preprocessed_blast_res/pires1200_bitscore_qcov_pid.tblastn6outmod"
  )


OUT_PATH <- here("data/processed/002.find_genes")

#Should I remove the excess of lone TrbH and traT on contigs?
REMOVE_LONE_TRBH_TRAT <- TRUE

#should I filter hits? preprocessed files already have been filtered to have only the top hits above 50% pid and 30% qcov
FILTER_HITS <- FALSE

#thresholds to use; applicable only if FILTER_HITS is TRUE
PID_THRESHOLD <- 50

QCOV_TRESHOLD <- 30


#################################################
# READ IN DATA


### Read in  blast results
pires1200_blast <- load_processed_blast6out(INPUT_BLAST_FILE)



#################################################
# FILTERING STEPS



if (FILTER_HITS) {
  pires1200_blast <-
    keep_good_hits(pires1200_blast,
                   pid_threshold = PID_THRESHOLD,
                   qcov_threshold = QCOV_TRESHOLD)
}




if (REMOVE_LONE_TRBH_TRAT) {
  #Summarise results to get one-hit contigs
  contig_summary_df <-
    summarise_at_contig_level(pires1200_blast)
  
  one_hit_contigs <- contig_summary_df %>%
    filter(n_hits_contig == 1) %>%
    pull(genome_contig)
  
  lone_trbh_contigs <-
    pires1200_blast %>% filter(genome_contig %in% one_hit_contigs) %>% filter(gene_name == "trbh") %>% pull(genome_contig)
  
  lone_trat_contigs <-
    pires1200_blast %>% filter(genome_contig %in% one_hit_contigs) %>% filter(gene_name == "trat") %>% pull(genome_contig)
  
  pires1200_blast <-
    pires1200_blast %>% filter(!genome_contig %in% lone_trbh_contigs,!genome_contig %in% lone_trat_contigs)
  
}


#################################################
# CREATE SUMMARY TABLES


contig_summary_df <-
  summarise_at_contig_level(pires1200_blast)


genome_summary_df <-
  summarise_at_genome_level(pires1200_blast)


#################################################
# SAVE FILTERED DATA AND SUMMARY TABLES


write.table(
  pires1200_blast,
  paste0(
    OUT_PATH,
    "/pires1200blast_pid_",
    PID_THRESHOLD,
    "_qcov",
    QCOV_TRESHOLD,
    "_rm1hits",
    REMOVE_LONE_TRBH_TRAT,
    ".tsv"
  ),
  quote = F,
  sep = '\t'
)

write.table(
  contig_summary_df,
  paste0(
    OUT_PATH,
    "/contig_summary_pires1200blast_",
    PID_THRESHOLD,
    "_qcov",
    QCOV_TRESHOLD,
    "_rm1hits_",
    REMOVE_LONE_TRBH_TRAT,
    ".tsv"
  ),
  sep = '\t',
  quote = F
)

write.table(
  genome_summary_df,
  paste0(
    OUT_PATH,
    "/genome_summary_pires1200blast_",
    PID_THRESHOLD,
    "_qcov",
    QCOV_TRESHOLD,
    "_rm1hits_",
    REMOVE_LONE_TRBH_TRAT,
    ".tsv"
  ),
  sep = '\t',
  quote = F
)


