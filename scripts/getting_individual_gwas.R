#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

condition_defaults <- list(
  control = "/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP/individual_lipids_final/control/raw_gwas/filtered_gwas",
  lowinput = "/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP/individual_lipids_final/lowinput/raw_gwas/filtered_gwas"
)

default_config <- list(
  condition = "control",
  gwas_dir = condition_defaults$control,
  gff = "/rsstu/users/r/rrellan/DOE_CAREER/SAP/ref/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3",
  outdir = ".",
  annotate_p_threshold = 0.05,
  gene_p_threshold = 1e-7,
  flank_bp = 25000L
)

usage_text <- paste(
  "Usage:",
  "  Rscript getting_individual_gwas.R --condition control --outdir /path/to/SAP-Lipidomics-Database",
  "",
  "Optional arguments:",
  "  --gwas-dir <path>              Override the default input folder for the selected condition",
  "  --gff <path>                   GFF3 gene annotation file",
  "  --annotate-p-threshold <num>   SNP p-value cutoff for gene-overlap annotation (default 0.05)",
  "  --gene-p-threshold <num>       Gene MinPValue cutoff kept in the final RDS (default 1e-7; this is -log10(p) >= 7)",
  "  --flank-bp <int>               Bases upstream/downstream of each SNP for gene overlap (default 25000)",
  "",
  "Examples:",
  "  Rscript getting_individual_gwas.R --condition control --outdir /Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database",
  "  Rscript getting_individual_gwas.R --condition lowinput --outdir /Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database",
  sep = "\n"
)

parse_args <- function(args) {
  parsed <- list()
  i <- 1

  while (i <= length(args)) {
    key <- args[[i]]

    if (key %in% c("-h", "--help")) {
      parsed$help <- TRUE
      i <- i + 1
      next
    }

    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }

    if (i == length(args)) {
      stop("Missing value for argument: ", key, call. = FALSE)
    }

    value <- args[[i + 1]]
    parsed[[sub("^--", "", key)]] <- value
    i <- i + 2
  }

  parsed
}

merge_config <- function(defaults, overrides) {
  out <- defaults
  for (nm in names(overrides)) {
    out[[nm]] <- overrides[[nm]]
  }
  out
}

normalize_condition <- function(x) {
  x <- tolower(trimws(x))
  if (!x %in% c("control", "lowinput")) {
    stop("--condition must be one of: control, lowinput", call. = FALSE)
  }
  x
}

get_trait_name <- function(path) {
  name <- basename(path)
  name <- sub("\\.assoc\\.txt$", "", name)
  name <- sub("\\.txt$", "", name)
  name <- sub("_mod_sub_Final_(control|lowinput)_all_lipids_BLUPs.*$", "", name)
  name
}

normalize_chr <- function(x) {
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  suppressWarnings(x_num <- as.integer(x))
  ifelse(is.na(x_num), x, x_num)
}

empty_annotation_df <- function() {
  data.frame(
    GeneID = character(),
    Chromosome = character(),
    `log(p)` = numeric(),
    stringsAsFactors = FALSE
  )
}

collapse_overlaps <- function(overlap_data) {
  if (nrow(overlap_data) == 0) {
    return(data.frame(
      GeneID = character(),
      Chromosome = character(),
      SNPs = character(),
      SNP_Positions = character(),
      Gene_Start = numeric(),
      Gene_End = numeric(),
      PValues = character(),
      Relation = character(),
      DistanceToGene = character(),
      MinPValue = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  overlap_data %>%
    group_by(GeneID) %>%
    summarise(
      Chromosome = as.character(Chromosome[1]),
      SNPs = paste(unique(SNP), collapse = ","),
      SNP_Positions = paste(unique(SNP_Position), collapse = ","),
      Gene_Start = Gene_Start[1],
      Gene_End = Gene_End[1],
      PValues = paste(unique(as.character(PValue)), collapse = ","),
      Relation = paste(unique(Relation), collapse = ","),
      DistanceToGene = paste(unique(as.character(DistanceToGene)), collapse = ","),
      MinPValue = min(as.numeric(PValue), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(MinPValue)
}

annotate_one_file <- function(file, genes_only, config) {
  trait_name <- get_trait_name(file)

  raw_df <- vroom::vroom(file, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::select(rs, chr, ps, p_wald) %>%
    transmute(
      SNP = ifelse(is.na(rs) | rs == "", paste0("SNP_", ps), as.character(rs)),
      Chromosome = normalize_chr(chr),
      Position = as.numeric(ps),
      PValue = as.numeric(p_wald)
    ) %>%
    filter(!is.na(Chromosome), !is.na(Position), !is.na(PValue))

  if (nrow(raw_df) == 0) {
    warning("Skipping empty or invalid GWAS file: ", basename(file), call. = FALSE)
    return(list(
      trait = trait_name,
      annotation_for_rds = empty_annotation_df(),
      n_snps = 0L,
      n_genes = 0L
    ))
  }

  annotate_df <- raw_df %>%
    filter(PValue <= config$annotate_p_threshold)

  if (nrow(annotate_df) == 0) {
    collapsed <- collapse_overlaps(data.frame())
  } else {
    snp_ranges <- GRanges(
      seqnames = Rle(as.character(annotate_df$Chromosome)),
      ranges = IRanges(start = annotate_df$Position, end = annotate_df$Position),
      marker = annotate_df$SNP,
      pvalue = annotate_df$PValue
    )

    expanded_snp_ranges <- resize(
      snp_ranges,
      width = as.integer(config$flank_bp) * 2L + 1L,
      fix = "center"
    )

    overlaps <- findOverlaps(genes_only, expanded_snp_ranges)

    if (length(overlaps) == 0) {
      collapsed <- collapse_overlaps(data.frame())
    } else {
      overlap_data <- data.frame(
        Chromosome = as.character(seqnames(expanded_snp_ranges[subjectHits(overlaps)])),
        GeneID = mcols(genes_only)[queryHits(overlaps), "ID"],
        SNP = mcols(expanded_snp_ranges)[subjectHits(overlaps), "marker"],
        SNP_Position = start(snp_ranges[subjectHits(overlaps)]),
        Gene_Start = start(genes_only[queryHits(overlaps)]),
        Gene_End = end(genes_only[queryHits(overlaps)]),
        PValue = mcols(snp_ranges)[subjectHits(overlaps), "pvalue"],
        stringsAsFactors = FALSE
      )

      overlap_data$Relation <- ifelse(
        overlap_data$Gene_Start <= overlap_data$SNP_Position & overlap_data$Gene_End >= overlap_data$SNP_Position,
        "within",
        ifelse(
          overlap_data$SNP_Position < overlap_data$Gene_Start,
          "upstream",
          "downstream"
        )
      )

      overlap_data$DistanceToGene <- ifelse(
        overlap_data$Relation == "within",
        "within",
        ifelse(
          overlap_data$Relation == "upstream",
          overlap_data$Gene_Start - overlap_data$SNP_Position,
          overlap_data$SNP_Position - overlap_data$Gene_End
        )
      )

      overlap_data$GeneID <- gsub("^gene:", "", overlap_data$GeneID)
      collapsed <- collapse_overlaps(overlap_data)
    }
  }

  annotation_for_rds <- collapsed %>%
    mutate(`log(p)` = -log10(MinPValue)) %>%
    filter(is.finite(`log(p)`), MinPValue <= config$gene_p_threshold) %>%
    transmute(
      GeneID = as.character(GeneID),
      Chromosome = as.character(Chromosome),
      `log(p)` = as.numeric(`log(p)`)
    )

  if (nrow(annotation_for_rds) == 0) {
    annotation_for_rds <- empty_annotation_df()
  }

  list(
    trait = trait_name,
    annotation_for_rds = annotation_for_rds,
    n_snps = nrow(raw_df),
    n_genes = nrow(annotation_for_rds)
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  cat(usage_text, sep = "\n")
  quit(save = "no", status = 0)
}

config <- merge_config(default_config, args)
config$condition <- normalize_condition(config$condition)

if (is.null(args$gwas_dir)) {
  config$gwas_dir <- condition_defaults[[config$condition]]
}

config$annotate_p_threshold <- as.numeric(config$annotate_p_threshold)
config$gene_p_threshold <- as.numeric(config$gene_p_threshold)
config$flank_bp <- as.integer(config$flank_bp)

if (!dir.exists(config$gwas_dir)) {
  stop("GWAS directory does not exist: ", config$gwas_dir, call. = FALSE)
}

if (!file.exists(config$gff)) {
  stop("GFF file does not exist: ", config$gff, call. = FALSE)
}

files <- list.files(config$gwas_dir, pattern = "\\.assoc\\.txt$|\\.txt$", full.names = TRUE)
if (length(files) == 0) {
  stop("No GWAS .txt files found in: ", config$gwas_dir, call. = FALSE)
}

ref_granges <- rtracklayer::import(config$gff)
genes_only <- ref_granges[mcols(ref_granges)$type == "gene"]

rds_dir <- file.path(config$outdir, "data", "gwas")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

results <- lapply(files, annotate_one_file, genes_only = genes_only, config = config)

annotation_list <- setNames(
  lapply(results, `[[`, "annotation_for_rds"),
  vapply(results, `[[`, character(1), "trait")
)

rds_path <- file.path(rds_dir, paste0("all_annotations_", config$condition, ".rds"))
saveRDS(annotation_list, rds_path)

unique_genes <- sort(unique(unlist(lapply(annotation_list, function(df) df$GeneID))))
write.table(
  unique_genes,
  file = file.path(rds_dir, paste0("unique_genes_", config$condition, ".txt")),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

manifest <- data.frame(
  trait = vapply(results, `[[`, character(1), "trait"),
  n_snps = vapply(results, `[[`, numeric(1), "n_snps"),
  n_genes_in_rds = vapply(results, `[[`, numeric(1), "n_genes"),
  stringsAsFactors = FALSE
)

write.csv(
  manifest,
  file = file.path(rds_dir, paste0("individual_gwas_manifest_", config$condition, ".csv")),
  row.names = FALSE
)
