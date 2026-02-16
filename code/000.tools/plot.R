## ---------------------------
##
## Script name: Plotting functions 
##
## Purpose of script:
## Standard plotting functions and themes 
##
##
## Author: Sneha Sundar
##
## Date Modified: 2023-08-08
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

##------------------------



## LIBRARIES
###############################

library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(paletteer)
library(svglite)
library(patchwork)
library("ggridges")
library(tidytext)
library(UpSetR)
library(ggalluvial)
library(ggparallel)

## Label_info

TRA_GENES_GROUP <- read.table(here("data/processed/002.find_genes/reference_tra_genes/TRA_GENES_LABELS.tsv"),header = T, quote = "", sep = '\t')
ESS_TRA <- TRA_GENES_GROUP %>% filter(Function != "nonessential", Function != "self-transfer prevention", TRA_GENES != "fino" ) %>% pull(TRA_GENES)

#write.table(TRA_GENES_GROUP, file = here("data/processed/002.find_genes/reference_tra_genes/TRA_GENES_LABELS.tsv"),quote = F,sep = "\t",col.names = T)
## THEME
###############################

#theme_font_axes <- theme(axis.text=element_text(size=24), 
                    #axis.title=element_text(size=22)) 


theme_font_axes <- theme(axis.text=element_text(size=26), 
                         axis.title=element_text(size=24)) 

theme_font_facet <-  theme(strip.text.x = element_text(size = 26, 
                                                       color = "black", face = "bold"))

theme_font_axes_small <- theme(axis.text=element_text(size=16), 
                         axis.title=element_text(size=20)) 

theme_font_small <-  theme(strip.text.x = element_text(size = 16, 
                                                       color = "black", face = "bold"))


theme_legend <- theme(legend.title = element_text(size = 24, face = "bold", color = "black"), # Change title font
legend.text = element_text(size = 22))

theme_legend_inside_right_top <- theme(legend.title = element_text(size = 16, face = "bold", color = "black"), # Change title font
                             legend.text = element_text(size = 14),
                             legend.position = "inside",
                             legend.position.inside = c(0.95,0.95),
                             legend.justification = c("right", "top"),
                             legend.box.just = "right",
                             legend.box.background = element_rect(),
                             legend.margin = margin(6, 6, 6, 6))

theme_legend_inside_right_bottom <- theme(legend.title = element_text(size = 16, face = "bold", color = "black"), # Change title font
      legend.text = element_text(size = 14),
      legend.position = c(0.95, 0.05),        # inside panel, bottom-right
      legend.justification = c("right", "bottom"),
      legend.background = element_rect(fill = alpha("white", 0.6), colour = NA),
      legend.key.size = unit(0.4, "cm"),
      legend.box.background = element_rect(),
      legend.margin = margin(6, 6, 6, 6))



themeBW <- theme_bw()  +
  theme_font_axes +
  theme_font_facet + 
  theme_legend 
  
#Manuscrupt themes

theme_legend_inside_right_bottom_manuscript <- theme(legend.title = element_text(size = 11), # Change title font
                                          legend.text = element_text(size = 10),
                                          legend.position = c(0.99, 0.01),        # inside panel, bottom-right
                                          legend.justification = c("right", "bottom"),
                                          legend.background = element_rect(),
                                          legend.key.size = unit(0.3, "cm"),
                                          legend.box.background = element_rect(colour = "grey"),
                                          legend.margin = margin(3, 3, 3, 3))

theme_legend_inside_left_bottom_manuscript <- theme(legend.title = element_text(size = 11), # Change title font
                                                     legend.text = element_text(size = 10),
                                                     legend.position = c(0.01, 0.01),        # inside panel, bottom-left
                                                     legend.justification = c("left", "bottom"),
                                                     legend.background = element_rect(),
                                                     legend.key.size = unit(0.3, "cm"),
                                                     legend.box.background = element_rect(colour = "grey"),
                                                     legend.margin = margin(3, 3, 3, 3))



theme_legend_inside_right_top_manuscript <- theme(legend.title = element_text(size = 11), # Change title font
                                       legend.text = element_text(size = 10),
                                       legend.position = "inside",
                                       legend.position.inside = c(0.99,0.99),
                                       legend.justification = c("right", "top"),
                                       legend.box.just = "right",
                                       legend.box.background = element_rect(colour = "grey"),
                                       legend.key.size = unit(0.3, "cm"),
                                       legend.margin = margin(3, 3, 3, 3))

theme_legend_inside_left_top_manuscript <- theme(legend.title = element_text(size = 11), # Change title font
                                                  legend.text = element_text(size = 10),
                                                  legend.position = "inside",
                                                  legend.position.inside = c(0.01,0.99),
                                                  legend.justification = c("left", "top"),
                                                  legend.box.just = "left",
                                                  legend.box.background = element_rect(colour = "grey"),
                                                  legend.key.size = unit(0.3, "cm"),
                                                  legend.margin = margin(3, 3, 3, 3))

themeBW_manuscript <- theme_bw(base_size = 10)

theme_manuscript <- themeBW_manuscript  +
  theme(
    plot.title      = element_text(size = 14, hjust = 0.5),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 11),
    legend.title    = element_text(size = 11,face="plain"),
    legend.text     = element_text(size = 10),
    strip.text      = element_text(size = 12,color = "black", face = "bold"),
    plot.caption    = element_text(size = 9, hjust = 0),
    plot.subtitle   = element_text(size = 14, hjust = 0.5),
    legend.box.background = element_rect(),
    strip.background = element_rect(fill = "white")
  )

themenull <- theme(axis.line = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   panel.background = element_blank(),
                   panel.border = element_blank(),
                   strip.background = element_blank(),
                   panel.grid = element_blank())

## COLOUR PALLETES
################################

CB15 <- c(
  '#4477AA',  # blue
  '#66CCEE',  # sky blue
  '#228833',  # green
  '#CCBB44',  # yellow
  '#EE6677',  # red
  '#AA3377',  # purple
  '#BBBBBB',  # gray
  '#000000',  # black
  '#661100',  # brown
  '#6699CC',  # steel blue
  '#888888',  # medium gray
  '#44AA99',  # teal
  '#117733',  # dark green
  '#332288',  # dark blue
  '#DDCC77'   # sand
)

## BINARY

binary_colours <- c("#F7FBFF", '#034e7b')
binary_colours2 <- c("FALSE" = "grey","TRUE" = "black")

## Sequential
#-----------------------


#https://colorbrewer2.org/?type=sequential&scheme=PuBu&n=8
PuBu<- c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb',
                  '#74a9cf','#3690c0','#0570b0','#034e7b')

## Diverging
#-----------------------

#https://colorbrewer2.org/?type=diverging&scheme=PiYG&n=8
PiYg <- c('#c51b7d','#de77ae','#f1b6da','#fde0ef',
                    '#e6f5d0','#b8e186','#7fbc41','#4d9221')


## Qualitative
#------------------------

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#colour brewer https://colorbrewer2.org/#type=qualitative&scheme=Set2&n=8                
#qualitative set2 ; not colour-blind friendly

Set2 <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3',
              '#a6d854','#ffd92f','#e5c494','#b3b3b3')

#from https://colorbrewer2.org/?type=qualitative&scheme=Dark2&n=8
#qualitative dark2 ; not colour-blind friendly 

#colours for the animals 
Dark2 <- c('#1b9e77','#d95f02','#7570b3','#e7298a',
                    '#66a61e','#e6ab02','#a6761d','#666666')

### Colours for categories ##

#animal_colours <- c(Bovine = "#31846B", Poultry = "#D55E00",Swine = "#918FCC")

animal_colours <- c(
  Bovine  = "#4b8c75",  # muted teal-green
  Poultry = "#c0703a",  # earthy rust orange
  Swine   = "#8d8bbd"   # soft lavender-gray
)

#source_niche_colours <- c("Livestock" = "#4477AA","Food"= "#DDCC77")
source_niche_colours <- c(
  "Livestock" = "#00678a",  # cool blue
  "Food"      = "#d5c36b"   # muted mustard
)

#resistance_presence <- c("Yes" = "#984464", "No" = "#666666")
resistance_presence <- c(
  "Yes" = "#984464",  # same plum tone as tra palette
  "No"  = "#8c8c8c"   # soft neutral gray
)

# phg_colours <- c(
#   'A'  = '#5aa392',
#   'B1' = '#e67a53',
#   'B2' = '#7f8db8',
#   'C'  = '#d178b2',
#   'D'  = '#92c14a',
#   'E'  = '#e6c82a',
#   'F'  = '#cdb383'
# )

phg_colours <- c(
  'A'  = '#e6a176',  
  'B1' = '#bb5a3c',  
  'B2' = '#7a8f55',  
  'C'  = '#b37b83',  
  'D'  = '#1f6f6f',  
  'E'  = '#e3a64d', 
  'F'  = '#8da0cb'   
)

n.distinct.genes.category_colours <-  c("[0,5)" = "#cdcdcd", "[5,23)" = "#9fc8c8","[23,35]" = "#1f6f6f")

functional_colours <- c("missing" = "#cdcdcd", "incomplete" = "#9fc8c8","complete" = "#1f6f6f")

coluzzi_colours <- c(
  "no.transfer.gene" = "#cdcdcd",  # neutral gray
  "mobless"          = "#d5c36b",  # matches Food (yellow-beige)
  "mob"              = "#e3a64d",  # warm amber
  "pdconj"           = "#c8dfb8",  # pale sage
  "conj"             = "#5b8a36"   # muted olive green
)

#colours for tra gene function


tra_gene_function_colours <- c("regulation/relaxosome" = "#e6a176",
  "pilus formation" = "#56641a", 
  "T4SS" = "#00678a",
  "mating pair stability" = "#c0affb",
  "nonessential" = "#cdcdcd",
  "self-transfer prevention" = "#984464")


dnds_colours <- c("pos" = "red", "neut" = "grey", "neg" = "black")


#discrete count colours
count_colours <- c('0'="#F7FBFF",'1'="#C6DBEF",'2'="#6BAED6",'3'="#2171B5",'4'="#08306B")

replicon_transfer_gene_colors <- c("transfer gene" = "#08306B","replicon" = "#e6a176")

has_transfer_operon_colours <- c("yes" = "#08306B","no" = "#cdcdcd")


#1f6f6f

## plots
################################
histogram <- function(data, variable, binwidth, facet, rows = 1, theme = themeBW){
  
  facet_layer <- if (!rlang::quo_is_null(rlang::enquo(facet))) facet_wrap(vars({{ facet }}),nrow = rows)
  
  p <- data %>% ggplot(aes({{variable}})) + geom_histogram(binwidth = binwidth) + facet_layer
  
  p <- data %>% ggplot(aes({{variable}})) + geom_histogram(binwidth = binwidth) + facet_layer
  return(p + theme)
}


