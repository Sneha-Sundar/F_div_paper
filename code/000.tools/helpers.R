## ---------------------------
##
## Script name: Helper functions
##
## 
## Author: Sneha Sundar
##

##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


## ------------------------
## PACKAGES
## ------------------------

library(here)
library(tidyverse)

library(easystats)
library(epitools)
library(emmeans)
library(dunn.test)
library(CooccurrenceAffinity) 


## ------------------------
## INPUT FUNCTIONS
## ------------------------


load_raw_blast6out <- function(path, genome, extension) {
  blastout_fp <- paste0(path, '/', genome, extension)
  if (file.exists(blastout_fp)) {
    return(read.table(
      blastout_fp,
      sep = '\t',
      col.names = c(
        "target_genome",
        "query",
        "target",
        "pid",
        "alnlen",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "evalue",
        "bitscore"
      ),
      quote = ""
    ))
  }
}


load_raw_bakta <- function(path, genome, extension) {
  fp <- paste0(path, '/', genome, extension)
  if (file.exists(fp)) {
    return(read.table(
      fp,
      sep = '\t',
      comment.char = "#",
      col.names = c(
        "contig",
        "type",
        "start",
        "stop",
        "strand",
        "locus_tag",
        "gene",
        "product",
        "dbxrefs"
      ),
      quote = ""
    ))
  }
}


load_tra_gene_metadata <- function(file_path) {
  query_metadata <-
    read.table(file_path,
               sep = '\t',
               quote = "",
               header = T)
  return(query_metadata)
}

#' Title
#'
#' @param path
#' @param genome_md_fp
#'
#' @return
#' @export
#'
#' @examples
load_processed_blast6out <-
  function(path, genome_meta_fp = "data/raw/pires_1200/pires_1200_ecoli_metadata.tsv") {
    blast_colnames <-
      read_lines(
        file = here(
          "data/processed/002.find_genes/pires1200/ecoli/preprocessed_blast_res/colnames.txt"
        )
      )
    
    tbl <-
      read.table(
        file = path,
        header = F,
        col.names = blast_colnames,
        quote = "",
        sep = '\t'
      )
    
    if (file.exists(here(genome_meta_fp))) {
      #add host genome metadata
      genome_meta <-
        read.table(
          here(genome_meta_fp),
          quote = "",
          sep = '\t',
          header = TRUE
        )
      
      tbl <- tbl %>%
        left_join(genome_meta %>%
                    select(
                      c(
                        "BioSample",
                        "Collection_Date",
                        "Latitude",
                        "Longitude",
                        "Generic_Host",
                        "Source_Niche",
                        "Source_Type",
                        "PhG",
                        "ST",
                        "CC",
                        "Number_Inc",
                        "Acq_Resist",
                        "Year"
                      )
                    ),
                  join_by(target_genome == BioSample))
      
      return(tbl)
      
    }
    return(tbl)
  }

load_processed_bakta <-
  function(path, genome_meta_fp = "data/raw/pires_1200/pires_1200_ecoli_metadata.tsv") {
    bakta_colnames <-
      read_lines(
        file = here(
          "data/processed/002.find_genes/pires1200/ecoli/preprocessed_bakta_res/colnames.txt"
        )
      )
    
    tbl <-
      read.table(
        file = path,
        header = F,
        col.names = bakta_colnames,
        quote = "",
        sep = '\t'
      )
    
    genome_contig <- paste0(tbl$genome, '_', tbl$contig)
    
    tbl <- cbind(genome_contig, tbl)
    
    if (file.exists(here(genome_meta_fp))) {
      #add host genome metadata
      genome_meta <-
        read.table(
          here(genome_meta_fp),
          quote = "",
          sep = '\t',
          header = TRUE
        )
      
      tbl <- tbl %>%
        left_join(genome_meta %>%
                    select(
                      c(
                        "BioSample",
                        "Collection_Date",
                        "Latitude",
                        "Longitude",
                        "Generic_Host",
                        "Source_Niche",
                        "Source_Type",
                        "PhG",
                        "ST",
                        "CC",
                        "Number_Inc",
                        "Acq_Resist",
                        "Year"
                      )
                    ),
                  join_by(genome == BioSample))
      
      
      
      return(tbl)
      
    }
    
    return(tbl)
  }

## ----------------------------
## DATA WRANGLING
## ----------------------------

#' Add tra gene metadata (key: query protein id) to each row of blast results file
#'
#' @return Data.Frame blast results with tra gene metadata
add_tra_gene_metadata <- function(blast_res_df, metadata_df) {
  #remove redundant metadata
  metadata_df <- metadata_df %>%
    select(
      -c(
        'Protein',
        'GenBank.Number',
        'gene',
        "protein_id",
        'ReplacedBy',
        'seq_len',
        'Status',
        'Title_old',
        'Organism_old',
        'protein_des'
      )
    )
  
  return(blast_res_df %>% left_join(metadata_df, by = c("query" = "new_protein_ids")) %>%
           mutate(genome_contig = paste0(target_genome, '_', target)))
}



#' Add query coverage to blast results.
#'
#' @param blast_res_df Blast results file with qend, qstart and query_protein_len columns
#'
#' @return Data.Frame Blast results file with qcov column
compute_qcov <- function(blast_res_df) {
  return(blast_res_df %>% mutate(qcov = round((((qend - qstart) + 1) / query_protein_len
  ) * 100, 3)))
}

#' Assign strand to blast results. If tstart > tend, strand is -.
#'
#' @param tstart
#' @param tend
#'
#' @return vector (new_tstart, new_tend, strand)
assign_strand <- function(tstart, tend) {
  if (is.na(tstart) | is.na(tend)) {
    return(c(tstart, tend, NA))
  }
  
  if (tstart > tend) {
    return(c(tend, tstart, '-'))
  } else{
    return(c(tstart, tend, '+'))
  }
}


## ----------------------------
## FILTERING 
## ----------------------------



#' Remove hits with low percentage id and query coverage. Thresholds set by user.
#'
#' @param data blast results file with pid and qcov columns.
#' @param pid_threshold Between 0 and 100. Default: 80
#' @param qcov_threshold Between 0 and 100. Default: 30
#'
#' @return filtered blast results
keep_good_hits <- function(data, pid_threshold, qcov_threshold) {
  res <- data %>% filter(pid > pid_threshold, qcov > qcov_threshold)
  print(
    paste(
      nrow(data) - nrow(res),
      " hits removed with pid <=",
      pid_threshold,
      " and qcov <= ",
      qcov_threshold
    )
  )
  return(res)
}


#' Pick the top n hits
#'
#' @param data blast results
#' @param on "genome" or "contig"
#' @param n Number of best hit (top hit, top 2 hits, top 3 hits)
#' @return Data.Frame containing top n hits (see column "rank_hit") of tra gene at contig level or genome level.
#' @note Based on highest pid, qcov, aln_len and lowest number of mismatches and gap openings. No tie handling.
#'
pick_top_n_hits <- function(data, on = "contig", n = 1) {
  print(paste0("Top ", n, " hits on ", on, " kept."))
  if (on == "contig") {
    tmp <- data %>%
      group_by(target, gene_name) %>%
      arrange(
        desc(bitscore),
        desc(qcov),
        desc(pid),
        desc(alnlen),
        mismatch,
        gapopen,
        .by_group = T
      ) %>%
      slice_head(n = n) %>%
      mutate(rank_hit = row_number())
    
    print(paste(nrow(data) - nrow(tmp), " hits removed."))
    
    return(tmp)
  } else{
    tmp <- data %>%
      group_by(target_genome, gene_name) %>%
      arrange(
        desc(bitscore),
        desc(qcov),
        desc(pid),
        desc(alnlen),
        mismatch,
        gapopen,
        .by_group = T
      ) %>%
      slice_head(n = n) %>%
      mutate(rank_hit = row_number())
    
    print(paste(nrow(data) - nrow(tmp), " hits removed."))
    
    return(tmp)
  }
  
}







## -------------------

## SUMMMARISE
## -------------------


#' Title
#'
#' @param blast_res 
#' @param contig_lengths_fp 
#' @param genome_meta_fp 
#'
#' @return
#' @export
#'
#' @examples
summarise_at_contig_level <-
  function(blast_res,
           contig_lengths_fp = here("data/processed/002.find_genes/pires1200/ecoli/contig_lengths/pires_1200.contig_len"),
           genome_meta_fp = here("data/raw/pires_1200/pires_1200_ecoli_metadata.tsv")) {
    
    
    contig_summary_df <-
      blast_res %>% group_by(target_genome, target) %>% summarise(
        n_hits_contig = n(),
        median_pid = median(pid),
        median_cov = median(qcov)
      )
    
    if (file.exists(here(contig_lengths_fp))) {
      contig_lens_df <-
        read.table(
          here(contig_lengths_fp),
          col.names = c('genome', 'contig', 'contig_seq_len'),
          quote = ""
        )
      
      #contig length info
      contig_summary_df <-
        contig_summary_df %>% left_join(contig_lens_df,
                                        join_by(target_genome == genome, target == contig))
    }
    
    if (file.exists(here(genome_meta_fp))) {
      #add host genome metadata
      genome_meta <-
        read.table(
          here(genome_meta_fp),
          quote = "",
          sep = '\t',
          header = TRUE
        )
      
      contig_summary_df <- contig_summary_df %>%
        left_join(genome_meta %>%
                    select(
                      c(
                        "BioSample",
                        "Collection_Date",
                        "Latitude",
                        "Longitude",
                        "Generic_Host",
                        "Source_Niche",
                        "Source_Type",
                        "PhG",
                        "ST",
                        "CC",
                        "Number_Inc",
                        "Acq_Resist",
                        "Year"
                      )
                    ),
                  join_by(target_genome == BioSample))
      
    }
    
    contig_summary_df <-
      contig_summary_df %>% mutate(genome_contig = paste0(target_genome, '_', target))
    
    return(contig_summary_df)
    
  }

summarise_at_genome_level <- function(blast_res,genomes_list_fp = here("data/raw/pires_1200/BioSample_ecoli.txt"),
                                      genome_meta_fp = here("data/raw/pires_1200/pires_1200_ecoli_metadata.tsv")){
  
  
  genome_summary_df <- pires1200_blast %>% group_by(target_genome) %>% summarise(n_contigs_genome = n_distinct(target),n_hits_genome = n(),median_pid_genome = median(pid), median_cov_genome = median(qcov))
  

  
  if (file.exists(here(genome_meta_fp))) {
    #add host genome metadata
    genome_meta <-
      read.table(
        here(genome_meta_fp),
        quote = "",
        sep = '\t',
        header = TRUE
      )
    
   genome_summary_df <- genome_meta %>%
      select(
        c(
          "BioSample",
          "Collection_Date",
          "Latitude",
          "Longitude",
          "Generic_Host",
          "Source_Niche",
          "Source_Type",
          "PhG",
          "ST",
          "CC",
          "Number_Inc",
          "Acq_Resist",
          "Year"
        )
      ) %>% left_join(genome_summary_df,join_by(BioSample == target_genome)) %>% 
     rename(target_genome = BioSample) %>% 
     mutate(n_contigs_genome = replace_na(n_contigs_genome,0),n_hits_genome = replace_na(n_hits_genome,0),median_pid_genome = replace_na(median_pid_genome,0),median_cov_genome = replace_na(median_cov_genome,0))
    

}

  return(genome_summary_df)
}


## -------------------

## COMPUTE
## -------------------




#' Create presence absence matrix of whether hits on a contig are present in a specified gene list.
#'
#' @param contig_id Contig id
#' @param blast_results blast results of tra operon reference genes vs genomes. Contig_id should be called "genome_contig" and gene name should be in column called "gene_name"
#' @param gene_list List of genes to check for presence or absence.
#' @return A binary vector of length(gene_list) indicating presence and absence of each gene.
#' @examples
#'
#'
is_gene_set_in_contig <-
  function(contig_id, blast_results, gene_list) {
    hits_in_contig <- blast_results %>%
      filter(genome_contig %in% contig_id) %>%
      distinct(gene_name) %>%
      pull(gene_name)
    
    bin <- gene_list %in% hits_in_contig
    names(bin) <- gene_list
    return(bin * 1)
  }

#' Create presence absence matrix of whether hits on a genome are present in a specified gene list.
#'
#' @param genome_id genome id
#' @param blast_results blast results of tra operon reference genes vs genomes (can be top n hits or raw results; it does not matter for this function) . Should contain blast6out columns and metadata of query and target. Contig_id should be called "target" and gene name should be in column called "gene_name"
#' @param gene_list List of genes to check for presence or absence.
#' @returns A binary vector of length(gene_list) indicating presence and absence of each gene.
#' @examples
#'
#'
is_gene_set_in_genome <-
  function(genome_id, blast_results, gene_list) {
    hits_in_genome <- blast_results %>%
      filter(target_genome %in% genome_id) %>%
      distinct(gene_name) %>%
      pull(gene_name)
    
    
    bin <- gene_list %in% hits_in_genome
    names(bin) <- gene_list
    return(bin * 1)
  }

count_tra_genes_in_genome <-
  function(genome_id, blast_results, gene_list) {
    hits_in_genome <- blast_results %>%
      filter(target_genome %in% genome_id) %>%
      pull(gene_name) %>% table()
    
  
    
    bin <- rep(0,36)
    names(bin) <- gene_list
    bin[names(hits_in_genome)] <- hits_in_genome
    
    return(bin)
    
  }

count_tra_genes_in_genome_bakta <-
  function(genome_id, bakta_results, gene_list) {
    hits_in_genome <- bakta_results %>%
      filter(genome %in% genome_id) %>%
      pull(gene) 
    
    hits_in_genome <-table(hits_in_genome[hits_in_genome %in% gene_list])
    
    bin <- rep(0,36)
    names(bin) <- gene_list
    bin[names(hits_in_genome)] <- hits_in_genome
    
    return(bin)
    
  }


count_tra_genes_in_contig_bakta <-
  function(genome_contig_id, bakta_results, gene_list) {
    hits_in_genome_contig <- bakta_results %>%
      filter(genome_contig %in% genome_contig_id) %>% 
      pull(gene) 
    
    hits_in_genome_contig <-table(hits_in_genome_contig[hits_in_genome_contig %in% gene_list])
    
    bin <- rep(0,36)
    names(bin) <- gene_list
    bin[names(hits_in_genome_contig)] <- hits_in_genome_contig
    
    return(bin)
    
  }

get_tra_gene_contigs <- function(blast_results) {
  #works only when a single genome is there in the results
  return(unique(blast_results$target))
}

get_bakta_contigs <- function(bakta_results, contig_list) {
  return(bakta_results %>% filter(contig %in% contig_list))
}


closest_distance_to_edge <- function(tstart, tend, contig_len) {
  if ((tstart - 1) < (contig_len - tend)) {
    return(tstart - 1)
  } else{
    return(contig_len - tend)
  }
}

compute_tra_gene_prevalence <- function(mat, denom) {
  if(nrow(mat) == 0){
    return(NA)
  }else{return(colSums(mat) / denom)}
}

compute_tra_gene_completeness <-
  function(mat, total_tra_genes = 36) {
    return(rowSums(mat) / total_tra_genes)
  }

#' Compute the intergenic distance of genes in a contig given sorted vectors of start and end coordinates.
#'
#' @return
#' @export
#'
#' @examples
compute_intergenic_distance <- function(start_vector, end_vector) {
  if (length(start_vector) != 1 &&
      length(start_vector) == length(end_vector)) {
    return(start_vector[2:length(start_vector)] - end_vector[1:length(end_vector) -
                                                               1])
  } else{
    return(NULL)
  }
  
}


pull_contig_values <- function(df, key, key_value, colname) {
  return(df %>% filter(.data[[key]] == key_value) %>% pull(.data[[colname]]))
  
}



compute_contig_intergenic_distances <-
  function(contig_id,
           df,
           start_col = "start",
           end_col = "stop") {
    df <- df %>% arrange(.data[[start_col]])
    starts <-
      pull_contig_values(df, "genome_contig", contig_id, start_col)
    ends <- pull_contig_values(df, "genome_contig", contig_id, end_col)
    compute_intergenic_distance(starts, ends)
  }

#' Given a presence absence pattern that is a character of 0s and 1s, convert it to a vector. 
#'
#' @return
#' @export
#'
#' @examples

presence_absence_pattern_to_vec <- function(pattern){
  
  as.numeric(strsplit(pattern, "")[[1]])
}


#' Given a presence absence pattern, return a list of missing genes . 
#'
#' @return
#' @export
#'
#' @examples


return_missing_genes <- function(pattern, gene_list = TRA_GENES_GROUP$TRA_GENES){
  
  vec <- presence_absence_pattern_to_vec(pattern)
  
  TRA_GENES_GROUP$TRA_GENES[!vec]
}


## ----------------------------
## PLOTTING 
## ----------------------------

#' Title
#'
#' @param blast_data
#' @param pid_threshold
#' @param qcov_threshold
#'
#' @return
#' @export
#'
#' @examples
plot_pid_vs_qcov <-
  function(blast_data,
           pid_threshold,
           qcov_threshold) {
    p <-
      ggplot() + geom_point(
        data = blast_data %>% filter(pid <= pid_threshold, qcov <= qcov_threshold),
        aes(x = pid, y = qcov, fill = target),
        shape = 21,
        color = "black",
        alpha = 0.4,
        size = 2,
        position = "jitter"
      ) + geom_point(
        data = blast_res %>% filter(pid > pid_threshold, qcov > qcov_threshold),
        aes(x = pid, y = qcov, fill = target),
        shape = 21,
        color = "black",
        alpha = 1,
        size = 2,
        position = "jitter"
      ) + xlab("Percentage identity") + ylab("Coverage") + theme_classic() + theme(axis.text =
                                                                                     element_text(size = 12),
                                                                                   axis.title = element_text(size = 14))
    
    return(p)
    
  }



has_tra_gene_in_contig_bakta <- function(contig_id, genome_bakta_df, gene_list = TRA_GENES){
  
  #check if any of the genes in the contig is a tra gene (based on annotation)
  return(any((genome_bakta_df %>% filter(contig == contig_id) %>% pull(gene)) %in% gene_list))
  
  
}


extract_from_label <- function(label, value_to_extract = "genome_contig"){
  
  
  values <- strsplit(label, "_")[[1]]
  
  if(value_to_extract == "genome"){
    return(values[1])
  }else if(value_to_extract == "contig"){
    paste(values[2],values[3],sep = "_")
  }else if(value_to_extract == "genome_contig"){
    paste(values[1],values[2],values[3],sep = "_")
  }else if(value_to_extract == "start"){
    as.numeric(values[4])
  }else if(value_to_extract == "end"){
    as.numeric(values[5])
  }else if(value_to_extract == "strand"){
    values[5]
  }
  

  
  
}

is_functional <- function(pres_abs_vec, req_genes, threshold){
  
  #check if at least `threshold` number of req_genes is present in a named pres_absence vector of genes
  #req genes must match the names in the pres_absence vector for this to work 
  
  return(sum(pres_abs_vec[req_genes] != 0) >= threshold)
  
}


create_prev_table <- function(group_cat = "Generic_Host", hit_type = "distinct_functional", min.hit, max.hit, genome_summary_table, PresAbsMat){
  
  
  prev_list = list()
  sample_size = numeric()
  
  if(group_cat %in% colnames(genome_summary_table)){
    
    for(group in  unique(genome_summary_table[[group_cat]])){
      
      if(hit_type == "distinct_functional"){
        
        group_genomes <- genome_summary_table %>% filter(.data[[group_cat]] == group, n.distinct.func.genes >= min.hit, n.distinct.func.genes <= max.hit) %>% pull(target_genome)
        
        group.prev <- compute_tra_gene_prevalence(PresAbsMat[group_genomes,],length(group_genomes))
        
        prev_list[[group]] <- group.prev
        
        sample_size <- c(sample_size, length(group_genomes))
      }
      
      if(hit_type == "functional"){
        
        group_genomes <- genome_summary_table %>% filter(.data[[group_cat]] == group, n.func.genes >= min.hit, n.func.genes <= max.hit) %>% pull(target_genome)
        
        group.prev <- compute_tra_gene_prevalence(PresAbsMat[group_genomes,],length(group_genomes))
        
        prev_list[[group]] <- group.prev
        
        sample_size <- c(sample_size, length(group_genomes))
      }
      
      
      if(hit_type == "all"){
        
        group_genomes <- genome_summary_table %>% filter(.data[[group_cat]] == group, bakta.n_hits_genome >= min.hit, bakta.n_hits_genome <= max.hit) %>% pull(target_genome)
        
        group.prev <- compute_tra_gene_prevalence(PresAbsMat[group_genomes,],length(group_genomes))
        
        prev_list[[group]] <- group.prev
        
        sample_size <- c(sample_size, length(group_genomes))
        
      }
      
      
      if(hit_type == "distinct"){
        
        group_genomes <- genome_summary_table %>% filter(.data[[group_cat]] == group, bakta.n_distinct_hits >= min.hit, bakta.n_distinct_hits <= max.hit) %>% pull(target_genome)
        
        group.prev <- compute_tra_gene_prevalence(PresAbsMat[group_genomes,],length(group_genomes))
        
        prev_list[[group]] <- group.prev
        
        sample_size <- c(sample_size, length(group_genomes))
      }
      
    }}else{
      if(hit_type == "distinct_functional"){
        
        group_genomes <- genome_summary_table %>% filter( n.distinct.func.genes >= min.hit, n.distinct.func.genes <= max.hit) %>% pull(target_genome)
        
        group.prev <- compute_tra_gene_prevalence(PresAbsMat[group_genomes,],length(group_genomes))
        
        prev_list[["all"]] <- group.prev
        
        sample_size <- c(sample_size, length(group_genomes))
      }
      
      if(hit_type == "functional"){
        
        group_genomes <- genome_summary_table %>% filter(n.func.genes >= min.hit, n.func.genes <= max.hit) %>% pull(target_genome)
        
        group.prev <- compute_tra_gene_prevalence(PresAbsMat[group_genomes,],length(group_genomes))
        
        prev_list[["all"]] <- group.prev
        
        sample_size <- c(sample_size, length(group_genomes))
      }
      
      
      if(hit_type == "all"){
        
        group_genomes <- genome_summary_table %>% filter( bakta.n_hits_genome >= min.hit, bakta.n_hits_genome <= max.hit) %>% pull(target_genome)
        
        group.prev <- compute_tra_gene_prevalence(PresAbsMat[group_genomes,],length(group_genomes))
        
        prev_list[["all"]] <- group.prev
        
        sample_size <- c(sample_size, length(group_genomes))
        
      }
      
      
      if(hit_type == "distinct"){
        
        group_genomes <- genome_summary_table %>% filter(bakta.n_distinct_hits >= min.hit, bakta.n_distinct_hits <= max.hit) %>% pull(target_genome)
        
        group.prev <- compute_tra_gene_prevalence(PresAbsMat[group_genomes,],length(group_genomes))
        
        prev_list[["all"]] <- group.prev
        
        sample_size <- c(sample_size, length(group_genomes))
      }
      
    }
  
  
  
  
  result <- list()
  result[["df"]] <- prev_list
  
  
  names(sample_size) <- unique(genome_summary_table[[group_cat]])
  result[["sample_size"]] <- sample_size
  return(result)
}



## -------------------

## SIMULATE
## -------------------


#' Sample n from a list of categories and output a presence absence vector of length given by the number of categories 
#'



sample_categories <- function(n, labels, re = F, weights){
  
  pres_abs_vec <- rep(0, length(labels))
  names(pres_abs_vec) <- labels
  categories <- sample(x = labels,size = n,replace = re ,prob = weights)
  
  pres_abs_vec[categories] <-1
  
  return(pres_abs_vec)
}


## -------------------

## Statistics
## -------------------



extract_lrt <- function(model, full_model) {
  a <- anova(full_model, model, test = "LRT")
  
  data.frame(
    df = a$Df[2],
    deviance = a$Deviance[2],
    residual_df = a$`Resid. Df`[2],
    residual_deviance = a$`Resid. Dev`[2],
    AIC = AIC(model),
    p_value = a$`Pr(>Chi)`[2]
  )
}
