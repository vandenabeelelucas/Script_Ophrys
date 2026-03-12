# Rscript: data_to_twisst.R
# Usage:
# Rscript data_to_twisst.R \
#   --input_pos file.tsv \
#   --tree_folder /path/to/trees \
#   --tree_pattern treefile \
#   --taxa H1,H2,H3,Outgroup \
#   --prune_tree TRUE/FALSE \
#   --chromo_size file.tsv
#   --alignment_size file.tsv

library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
library(optparse)

option_list <- list(
  make_option(c("--ali_position"), type="character", default=NULL,
              metavar="character"),
  make_option(c("--tree_folder"), type="character", default=NULL,
              metavar="character"),
  make_option(c("--tree_pattern"), type="character", default=NULL,
              metavar="character"),
  make_option(c("--taxa"), type="character", default=NULL,
              metavar="character"),
  make_option(c("--prune_tree"), type="logical", default=TRUE,
              metavar="logical"),
  make_option(c("--chromosome_size"), type="character", default=NULL,
              metavar="character"),
  make_option(c("--alignment_size"), type="character", default=NULL,
              metavar="character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Input files and parameters
input_pos            <- opt$ali_position
tree_folder          <- opt$tree_folder
tree_pattern         <- opt$tree_pattern
taxa                 <- strsplit(opt$taxa, ",")[[1]]
prune_tree           <- opt$prune_tree
chromo_size_file     <- opt$chromosome_size
alignment_size_file  <- opt$alignment_size

### READ ALIGNMENT POSITIONS ON REFERENCE GENOME
  
# Read input position file (tab-separated table) containing
# genomic coordinates for each alignment. The file must include:
#   - "ali"   : alignment identifier (alignment name)
#   - "chrom" : chromosome name
#   - "start" : start position of the alignment on the reference genome
#   - "end"   : end position of the alignment on the reference genome

pos <- read.table(input_pos,
                 header = TRUE,
                 sep = "\t",
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

# Keep only relevant columns
pos <- pos[, c("ali", "chrom", "start", "end")]

# Convert genomic coordinates to numeric
pos$start <- as.numeric(pos$start)
pos$end   <- as.numeric(pos$end)

# Ensure start < end
swap_idx <- which(pos$start > pos$end)
if (length(swap_idx) > 0) {
  temp <- pos$start[swap_idx]
  pos$start[swap_idx] <- pos$end[swap_idx]
  pos$end[swap_idx] <- temp
}

#Split alignments by chromosome
pos_by_chrom <- split(pos, pos$chrom)

#Sort each chromosome by start position
pos_by_chrom <- lapply(pos_by_chrom, arrange, start)

#Ensure no overlaps between alignments
resolve_overlaps <- function(df) {
  for (i in 2:nrow(df)) {
    if (df$end[i - 1] >= df$start[i]) {
      df$end[i - 1] <- df$start[i] - 1
      if (df$end[i - 1] < df$start[i - 1]) {
        df$end[i - 1] <- df$start[i - 1]
      }
    }
  }
  return(df)
}
pos_by_chrom <- lapply(pos_by_chrom, resolve_overlaps)

### TOPOLOGY ASSIGNMENT FOR EACH ALIGNMENTS

#Load trees
listinitial <- list.files(path = tree_folder, pattern = paste0("\\.", tree_pattern, "$"), full.names = TRUE)

list_tree <- lapply(listinitial, read.tree)

#Function for topology assignment
Topo_assignation <- function(treeset, list_files, taxa, prune_tree) {
  
  result <- data.frame(
    ali = basename(list_files),
    topo1 = integer(length(treeset)),
    topo2 = integer(length(treeset)),
    topo3 = integer(length(treeset)),
    notopo = integer(length(treeset)),
    stringsAsFactors = FALSE
  )
  
  # Prune trees if requested
  if (prune_tree) {
    pruned_trees <- lapply(treeset, function(tr) keep.tip(tr, taxa))
  } else {
    pruned_trees <- treeset
  }
   class(pruned_trees) <- "multiPhylo"
   
   # Root using outgroup 
  pruned_trees <- root(pruned_trees, outgroup = taxa[4], resolve.root = TRUE)
  
  # Assign topology
  for (i in seq_along(pruned_trees)) {
    tree_to_test <- pruned_trees[[i]]
    if (is.monophyletic(tree_to_test, taxa[1:2])) {
      result[i, 2:5] <- c(1,0,0,0)
    } else if (is.monophyletic(tree_to_test, taxa[2:3])) {
      result[i, 2:5] <- c(0,1,0,0)
    } else {
      result[i, 2:5] <- c(0,0,1,0)
    }
  }
  
  return(result)
}

#Launch Topology assignation
topo <- Topo_assignation(list_tree, listinitial, taxa, prune_tree = prune_tree)

#Merge topology with genomic positions
topo_pos <- topo %>%
  left_join(pos[, c("ali", "chrom", "start", "end")],
            by = "ali")
topo_by_chrom <- split(topo_pos, topo_pos$chrom)
topo_by_chrom <- lapply(topo_by_chrom, arrange, start)

### BUILD PSEUDO-CHROMOSOMES

#Load alignment size
# The file must include:
#   - "ali"   : alignment identifier (alignment name), same as input_pos
#   - "size" : alignment size

size <- read.table(alignment_size_file,
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  check.names = FALSE)

#"Pseudo-chromosomes" generation
pos_by_pchrom <- lapply(pos_by_chrom, function(df) {
  
  df <- left_join(df, size[, c("ali", "size")], by = "ali") |>
    select(-any_of(c("start","end"))) |>
    mutate(size = as.numeric(size))
  
  df$start <- 1 + cumsum(c(0, head(df$size + 1, -1)))
  df$end   <- df$start + df$size
  
  df
})

### ADD NO_DATA SEGMENTS ALONG CHROMOSOMES

#Load chromosome full size
# The file must include:
#   - "chrom" : chromosome name
#   - "chrom_size" : chromosome size

chrom_sizes <- read.table(chromo_size_file, 
                          header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)

chrom_sizes$chrom <- trimws(as.character(chrom_sizes$chrom))
chrom_sizes$chrom_size  <- as.numeric(chrom_sizes$chrom_size)
names(pos_by_chrom) <- trimws(names(pos_by_chrom))

# Function to add no_data segments
add_no_data_topo <- function(df, chrom_size) {
  new_df <- data.frame(
    ali   = character(),
    chrom = character(),
    start = integer(),
    end   = integer(),
    topo1 = integer(),
    topo2 = integer(),
    topo3 = integer(),
    notopo = integer(),
    stringsAsFactors = FALSE
  )
  
  # Beginning gap
  if(df$start[1] > 1) {
    new_df <- rbind(new_df, data.frame(
      ali   = "No_data",
      chrom = df$chrom[1],
      start = 1,
      end   = df$start[1] - 1,
      topo1 = 0, topo2 = 0, topo3 = 0, notopo = 1,
      stringsAsFactors = FALSE
    ))
  }
  
  # Internal gaps
  for(i in seq_len(nrow(df))) {
    new_df <- rbind(new_df, df[i, ])
    
    if(i < nrow(df)) {
      gap_start <- df$end[i] + 1
      gap_end   <- df$start[i+1] - 1
      if(gap_start <= gap_end) {
        new_df <- rbind(new_df, data.frame(
          ali   = "No_data",
          chrom = df$chrom[i],
          start = gap_start,
          end   = gap_end,
          topo1 = 0, topo2 = 0, topo3 = 0, notopo = 1,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Final gap
  if(df$end[nrow(df)] < chrom_size) {
    new_df <- rbind(new_df, data.frame(
      ali   = "No_data",
      chrom = df$chrom[1],
      start = df$end[nrow(df)] + 1,
      end   = chrom_size,
      topo1 = 0, topo2 = 0, topo3 = 0, notopo = 1,
      stringsAsFactors = FALSE
    ))
  }
  
  return(new_df)
}

# Launch add_no_data_topo function
topo_by_chrom_nd <- list()
for(chr in chrom_sizes$chrom) {
  chrom_size <- chrom_sizes$chrom_size[chrom_sizes$chrom == chr]
  topo_by_chrom_nd[[chr]] <- add_no_data_topo(topo_by_chrom[[chr]], chrom_size)
}

###EXPORT CHROMOSOMES DATA

#Create directory
outdir <- "twisst_files"
dir.create(outdir, showWarnings = FALSE)

#Define topology header
header_lines <- c(
  "#topo1 ((H1,H2),H3);",
  "#topo2 ((H2,H3),H1);",
  "#topo3 ((other,other),other);",
  "topo1\ttopo2\ttopo3\tnotopo"
)


#Export one topology file with all chromosomes
all_df <- do.call(rbind, topo_by_chrom_nd)[ , c("topo1","topo2","topo3","notopo")]
outfile_all <- file.path(outdir, "Chrom_topologies.txt")
writeLines(header_lines, outfile_all)
write.table(all_df,
            file = outfile_all,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE)


#Define position header
header_line <- "chrom\tstart\tend"

#Export file
all_df <- do.call(rbind, topo_by_chrom_nd)[ , c("chrom","start","end")]
outfile_all <- file.path(outdir, "Chrom_positions.txt")
writeLines(header_line, outfile_all)
write.table(all_df,
            file = outfile_all,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE)

###CREATE AND EXPORT CHROMOSOMES DATA WITH ONLY GF

#keep only information for GF topologies
topo_by_chrom_nd_GF <- lapply(topo_by_chrom_nd, function(df) {
  df_mod <- df
  idx <- df_mod$topo1 == 1 | df_mod$topo3 == 1
  df_mod$topo1[idx] <- 0
  df_mod$topo3[idx] <- 0
  df_mod$notopo[idx] <- 1
  return(df_mod)
})

#Export file
all_df <- do.call(rbind, topo_by_chrom_nd_GF)[ , c("topo1","topo2","topo3","notopo")]
outfile_all <- file.path(outdir, "Chrom_GF_topologies.txt")
writeLines(header_lines, outfile_all)
write.table(all_df,
            file = outfile_all,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE)

###EXPORT PSEUDO-CHROMOSOMES DATA

#Create directory
pseudo_outdir <- "twisst_pseudo_files"
dir.create(pseudo_outdir, showWarnings = FALSE)

#Define topology header
header_lines <- c(
  "#topo1 ((H1,H2),H3);",
  "#topo2 ((H2,H3),H1);",
  "#topo3 ((other,other),other);",
  "topo1\ttopo2\ttopo3\tnotopo"
)

#Export file
all_df <- do.call(rbind, topo_by_chrom)[ , c("topo1","topo2","topo3","notopo")]
outfile_all <- file.path(pseudo_outdir, "all_pseudo_chromosomes_topologies.txt")
writeLines(header_lines, outfile_all)
write.table(all_df,
            file = outfile_all,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE)


#Define position header
header_line <- "chrom\tstart\tend"

#Export one topology file with all pseudo-chromosome
all_df <- do.call(rbind, pos_by_pchrom)[ , c("chrom","start","end")]
outfile_all <- file.path(pseudo_outdir, "all_pseudo_chromosomes_positions.txt")
writeLines(header_line, outfile_all)
write.table(all_df,
            file = outfile_all,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE)

###LAUNCH TWISST PLOT:
# Load TWISST plotting functions
source("plot_twisst_modify.R")

#Load data

# Chromosomes
weights_chr <- file.path(outdir, "Chrom_topologies.txt")
weights_chr_GF <- file.path(outdir, "Chrom_GF_topologies.txt")
positions_chr <- file.path(outdir, "Chrom_positions.txt")

# Pseudo-chromosomes
weights_pchr <- file.path(pseudo_outdir, "all_pseudo_chromosomes_topologies.txt")
positions_pchr <- file.path(pseudo_outdir, "all_pseudo_chromosomes_positions.txt")

#Import data

# Chromosome data
twisst_chr <- import.twisst(
  weights_files = weights_chr,
  window_data_files = positions_chr
)

# Chromosome data with GF only
twisst_chr_GF <- import.twisst(
  weights_files = weights_chr_GF,
  window_data_files = positions_chr
)

# Pseudo-chromosome data
twisst_pchr <- import.twisst(
  weights_files = weights_pchr,
  window_data_files = positions_pchr
)

# Plot chromosomes (mode = 2)

jpeg(
  filename = "twisst_chromosomes.jpg",
  width = 7200,
  height = 4800,
  res = 200
)

plot.twisst(
  twisst_chr,
  mode = 2,
  show_topos = FALSE,
  include_region_names = TRUE
)

dev.off()

# Plot chromosomes GF only (mode = 2)

jpeg(
  filename = "twisst_chromosomes_GF.jpg",
  width = 7200,
  height = 4800,
  res = 200
)

plot.twisst(
  twisst_chr_GF,
  mode = 2,
  show_topos = FALSE,
  include_region_names = TRUE
)

dev.off()

# Plot pseudo_chromosomes (mode = 3)

jpeg(
  filename = "twisst_pseudochromosomes.jpg",
  width = 7200,
  height = 4800,
  res = 200
)

plot.twisst(
  twisst_pchr,
  mode = 3,
  show_topos = FALSE,
  include_region_names = TRUE
)

dev.off()

### PLOT TOPOLOGY PROPORTION:

#Calculate topology proportion in each chromosomes 
topo_proportion <- topo_pos %>%
  summarise(
    topo_1  = sum(topo1),
    topo_2  = sum(topo2),
    topo_3  = sum(topo3),
    .by = chrom
  ) %>%
  pivot_longer(cols = topo_1:topo_3,
               names_to = "topology",
               values_to = "count") %>%
  group_by(chrom) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

#Plot topology proportion
jpeg("Proportion_plot.jpg", width = 8*300, height = 6*300, res = 300)
ggplot(topo_proportion, aes(x = chrom, y = freq, fill = topology)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count),
            position = position_stack(vjust = 0.5),
            color = "white", size = 3) +
  scale_fill_manual(values = c(
    "topo_1" = "#0075DC",
    "topo_2"   = "#2BCE48",
    "topo_3" = "#FFA405"
  )) +
  labs(x = "", y = "") +   
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")
dev.off()


### PLOT GENOME COVERAGE:

# Prepare plot
nchr <- length(chrom_sizes$chrom)
ncol <- 3
nrow <- ceiling(nchr / ncol)
jpeg("Coverage.jpg",
     width = 7200, height = 9600, res = 200)
par(mfrow = c(nrow, ncol), mar = c(5,5,3,2), 
    cex.main = 2.5, cex.lab = 2, cex.axis = 1.8)

# Plot coverage for 1Mb window
for(chr in chrom_sizes$chrom) {
  genes_chr <- pos_by_chrom[[chr]]
  max_chr <- chrom_sizes$chrom_size[chrom_sizes$chrom == chr]
  
  hist(genes_chr$start,
       breaks = seq(0, max_chr + 1e6, by = 1e6),
       right = TRUE,
       main = chr,
       xlab = "",
       ylab = "",
       col = "steelblue",
       border = "white")
}

dev.off()
