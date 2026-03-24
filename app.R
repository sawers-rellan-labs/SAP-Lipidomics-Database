 library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(dplyr)
library(vroom)
library(plotly)
library(Rtsne)
library(umap)
library(heatmaply)
library(RColorBrewer)
library(tibble)
library(VennDiagram)
library(grid)
library(ggnewscale)
library(stringr)
library(shinycssloaders)  # Loading spinners
library(visNetwork)       # Network visualization
library(igraph)           # Graph/network analysis
library(viridisLite)      # Colorblind-friendly palettes

app_title <- "SAP Lipidomics Database"
data_type_choices <- c(
  "Raw Signals" = "raw",
  "%TIC (Normalized)" = "tic",
  "CLT/CLR (Compositional)" = "clr"
)

resolve_existing_path <- function(candidates) {
  matches <- candidates[file.exists(candidates)]
  if (length(matches) == 0) {
    return(NA_character_)
  }
  matches[[1]]
}

dependency_specs <- list(
  lowinput = list(
    label = "Low-input lipid matrix",
    required = TRUE,
    candidates = c(
      file.path("data", "Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv"),
      file.path("data", "Lowinput_all_lipids_final.csv"),
      file.path("data", "lowinput_lipids.csv"),
      file.path("data", "raw", "Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv")
    )
  ),
  control = list(
    label = "Control lipid matrix",
    required = TRUE,
    candidates = c(
      file.path("data", "Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv"),
      file.path("data", "Control_all_lipids_final.csv"),
      file.path("data", "control_lipids.csv"),
      file.path("data", "raw", "Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv")
    )
  ),
  lipid_class = list(
    label = "Lipid class mapping table",
    required = TRUE,
    candidates = c(
      file.path("data", "final_lipid_classes.csv"),
      file.path("data", "lipid_class.csv"),
      file.path("data", "metadata", "final_lipid_classes.csv")
    )
  ),
  gene_annotation = list(
    label = "Gene annotation table",
    required = TRUE,
    candidates = c(
      file.path("data", "annotation", "gene_annotation.txt"),
      file.path("data", "gene_annotation.txt"),
      file.path("data", "metadata", "gene_annotation.txt")
    )
  ),
  ann_lowinput = list(
    label = "GWAS annotations for lowinput lipids (-log10(p) >= 7)",
    required = FALSE,
    candidates = c(
      file.path("data", "GWAS_RDS", "individual", "all_annotations_lowinput_individual_plog107.rds"),
      file.path("data", "all_annotations_lowinput.rds"),
      file.path("data", "gwas", "all_annotations_lowinput.rds"),
      file.path("data", "GWAS_RDS", "individual", "all_annotations_lowinput_individual.rds")
    )
  ),
  ann_control = list(
    label = "GWAS annotations for control lipids (-log10(p) >= 7)",
    required = FALSE,
    candidates = c(
      file.path("data", "GWAS_RDS", "individual", "all_annotations_control_individual_plog107.rds"),
      file.path("data", "all_annotations_control.rds"),
      file.path("data", "gwas", "all_annotations_control.rds"),
      file.path("data", "GWAS_RDS", "individual", "all_annotations_control_individual.rds")
    )
  ),
  ann_lowinput_sum_ratio = list(
    label = "GWAS annotations for lowinput sum-ratio traits (-log10(p) >= 7)",
    required = FALSE,
    candidates = c(
      file.path("data", "GWAS_RDS", "sum_ratio", "all_annotations_lowinput_sum_ratio_plog107.rds"),
      file.path("data", "GWAS_RDS", "sum_ratio", "all_annotations_lowinput_sum_ratio.rds")
    )
  ),
  ann_control_sum_ratio = list(
    label = "GWAS annotations for control sum-ratio traits (-log10(p) >= 7)",
    required = FALSE,
    candidates = c(
      file.path("data", "GWAS_RDS", "sum_ratio", "all_annotations_control_sum_ratio_plog107.rds"),
      file.path("data", "GWAS_RDS", "sum_ratio", "all_annotations_control_sum_ratio.rds")
    )
  ),
  manhattan_lowinput = list(
    label = "Manhattan plot images for lowinput lipids",
    required = FALSE,
    candidates = c(file.path("www", "manhattans", "lowinput"))
  ),
  manhattan_control = list(
    label = "Manhattan plot images for control lipids",
    required = FALSE,
    candidates = c(file.path("www", "manhattans", "control"))
  )
)

resolved_paths <- lapply(dependency_specs, function(spec) {
  resolve_existing_path(spec$candidates)
})

dependency_status <- tibble::tibble(
  key = names(dependency_specs),
  label = vapply(dependency_specs, function(spec) spec$label, character(1)),
  required = vapply(dependency_specs, function(spec) spec$required, logical(1)),
  candidates = vapply(
    dependency_specs,
    function(spec) paste(spec$candidates, collapse = " or "),
    character(1)
  ),
  path = unname(unlist(resolved_paths, use.names = FALSE))
)

missing_required_files <- dependency_status %>%
  filter(required, is.na(path))

missing_optional_files <- dependency_status %>%
  filter(!required, is.na(path))

app_has_required_data <- nrow(missing_required_files) == 0

empty_lipid_dataset <- tibble::tibble(Sample = character())
empty_lipid_class_info <- tibble::tibble(Lipids = character(), Class = character())
empty_gene_annotation <- tibble::tibble(GeneID = character())
empty_gwas_lipid_map <- tibble::tibble(
  Original_Lipid = character(),
  Class = character(),
  Subclass = character()
)

read_lipid_matrix <- function(path) {
  if (is.na(path)) {
    return(empty_lipid_dataset)
  }

  df <- vroom::vroom(path, show_col_types = FALSE)

  if (ncol(df) >= 4) {
    df <- df %>% dplyr::select(-c(2:4))
  }

  df
}

read_lipid_class_info <- function(path) {
  if (is.na(path)) {
    return(empty_lipid_class_info)
  }

  vroom::vroom(path, show_col_types = FALSE) %>%
    filter(!is.na(Class))
}

read_gene_annotation <- function(path) {
  if (is.na(path)) {
    return(empty_gene_annotation)
  }

  vroom::vroom(path, delim = "\t", show_col_types = FALSE)
}

read_gwas_lipid_map <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(empty_gwas_lipid_map)
  }

  df <- vroom::vroom(path, show_col_types = FALSE)
  if (nrow(df) == 0) {
    return(empty_gwas_lipid_map)
  }

  lipid_col <- dplyr::case_when(
    "Original_Lipid" %in% names(df) ~ "Original_Lipid",
    "Lipids" %in% names(df) ~ "Lipids",
    TRUE ~ NA_character_
  )
  class_col <- dplyr::case_when(
    "Class" %in% names(df) ~ "Class",
    "Mapped_Nar Class" %in% names(df) ~ "Mapped_Nar Class",
    "Mapped_Nar_Class" %in% names(df) ~ "Mapped_Nar_Class",
    "Primary_Class" %in% names(df) ~ "Primary_Class",
    TRUE ~ NA_character_
  )
  subclass_col <- dplyr::case_when(
    "Subclass" %in% names(df) ~ "Subclass",
    "Secondary_Class" %in% names(df) ~ "Secondary_Class",
    "SubClass" %in% names(df) ~ "SubClass",
    TRUE ~ NA_character_
  )

  if (is.na(lipid_col) || is.na(class_col)) {
    return(empty_gwas_lipid_map)
  }

  out <- df %>%
    transmute(
      Original_Lipid = as.character(.data[[lipid_col]]),
      Class = as.character(.data[[class_col]]),
      Subclass = if (!is.na(subclass_col)) as.character(.data[[subclass_col]]) else ""
    ) %>%
    mutate(
      Original_Lipid = trimws(Original_Lipid),
      Class = ifelse(is.na(Class), "", trimws(Class)),
      Subclass = ifelse(is.na(Subclass), "", trimws(Subclass))
    ) %>%
    filter(Original_Lipid != "") %>%
    distinct(Original_Lipid, .keep_all = TRUE)

  if (nrow(out) == 0) {
    return(empty_gwas_lipid_map)
  }

  out
}

read_annotation_bundle <- function(path) {
  if (is.na(path)) {
    return(list())
  }

  obj <- readRDS(path)
  if (!is.list(obj)) {
    return(list())
  }
  obj
}

merge_annotation_lists <- function(...) {
  annotation_lists <- list(...)
  annotation_lists <- annotation_lists[vapply(annotation_lists, is.list, logical(1))]
  annotation_lists <- annotation_lists[vapply(annotation_lists, length, integer(1)) > 0]

  if (length(annotation_lists) == 0) {
    return(list())
  }

  merged <- do.call(c, annotation_lists)
  merged[!duplicated(names(merged))]
}

merge_annotation_lists_with_source <- function(individual = list(), sum_ratio = list(), other_label = "other") {
  individual <- if (is.list(individual)) individual else list()
  sum_ratio <- if (is.list(sum_ratio)) sum_ratio else list()

  combined <- c(individual, sum_ratio)
  if (length(combined) == 0) {
    return(list(data = list(), source = setNames(character(0), character(0))))
  }

  source_map <- c(
    setNames(rep(other_label, length(individual)), names(individual)),
    setNames(rep("sum_ratio", length(sum_ratio)), names(sum_ratio))
  )

  keep <- !duplicated(names(combined))
  list(
    data = combined[keep],
    source = source_map[keep]
  )
}

build_setup_banner <- function() {
  if (app_has_required_data && nrow(missing_optional_files) == 0) {
    return(NULL)
  }

  render_dependency_list <- function(df) {
    tags$ul(
      lapply(seq_len(nrow(df)), function(i) {
        tags$li(
          tags$strong(df$label[[i]]),
          tags$br(),
          tags$code(df$candidates[[i]])
        )
      })
    )
  }

  tags$div(
    class = "setup-banner",
    tags$h4("Shiny app setup"),
    tags$p("Drop the files below into this repository and restart the app to unlock the full interface."),
    if (nrow(missing_required_files) > 0) {
      tags$div(
        tags$strong("Required files still missing:"),
        render_dependency_list(missing_required_files)
      )
    },
    if (nrow(missing_optional_files) > 0) {
      tags$div(
        tags$strong("Optional assets not found yet:"),
        tags$p("The app will run without these, but GWAS Manhattan images will show a placeholder message."),
        render_dependency_list(missing_optional_files)
      )
    }
  )
}

lowinput <- read_lipid_matrix(resolved_paths$lowinput)
control <- read_lipid_matrix(resolved_paths$control)

dataset_list <- list(
  lowinput = lowinput,
  control = control
)

lipid_class_info <- read_lipid_class_info(resolved_paths$lipid_class)
subsubclass_col <- dplyr::case_when(
  "Sub-subclass" %in% names(lipid_class_info) ~ "Sub-subclass",
  "Sub_subclass" %in% names(lipid_class_info) ~ "Sub_subclass",
  TRUE ~ NA_character_
)
pca_class_choices <- sort(unique(stats::na.omit(as.character(lipid_class_info$Class))))
if (length(pca_class_choices) == 0) {
  pca_class_choices <- character(0)
}
gwas_lipid_class_map <- list(
  control = read_gwas_lipid_map(file.path("data", "metadata", "lipid_with_classes_control.csv")),
  lowinput = read_gwas_lipid_map(file.path("data", "metadata", "lipid_with_classes_lowinput.csv"))
)

#---------------------------
# TIC Normalization Function (%TIC per sample)
#---------------------------
tic_normalize <- function(df, id_col = 1) {
  if (is.null(df) || ncol(df) <= 1) {
    return(df)
  }

  id_name <- if (is.numeric(id_col)) names(df)[id_col] else id_col

  if (length(id_name) == 0 || is.na(id_name)) {
    return(df)
  }

  # Coerce all lipid columns to numeric (handles character-encoded numbers)
  X <- df %>%
    dplyr::select(-dplyr::all_of(id_name)) %>%
    as.data.frame(stringsAsFactors = FALSE, check.names = FALSE)
  X <- as.data.frame(lapply(X, function(col) {
    if (is.numeric(col)) {
      return(as.numeric(col))
    }
    col_chr <- gsub(",", "", trimws(as.character(col)))
    suppressWarnings(as.numeric(col_chr))
  }), stringsAsFactors = FALSE, check.names = FALSE)

  if (ncol(X) == 0) {
    return(df)
  }

  Xmat <- as.matrix(X)

  # Calculate row sums (TIC)

  rs <- rowSums(Xmat, na.rm = TRUE)
  rs[rs == 0 | is.na(rs)] <- NA_real_

  # Convert to percentage
  Xpct <- sweep(Xmat, 1, rs, FUN = "/") * 100
  Xpct[is.na(Xpct)] <- 0

  # Combine with ID column
  out <- df %>% dplyr::select(dplyr::all_of(id_name))
  bind_cols(out, as.data.frame(Xpct, check.names = FALSE))
}

# Create TIC-normalized versions of datasets
dataset_list_tic <- list(
  lowinput = tic_normalize(dataset_list$lowinput),
  control = tic_normalize(dataset_list$control)
)

#---------------------------
# Plot Theme (from paper script)
#---------------------------
plot_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

#---------------------------
# Class colors (from paper script - exact colors)
#---------------------------
class_colors <- c(
  PC   = "#00441B",
  PA   = "#1B7837",
  PE   = "#41AB5D",
  PG   = "#78C679",
  PS   = "#C2E699",
  DG   = "#54278F",
  DGDG = "#F768A1",
  MG   = "#8941ED",
  MGDG = "#FBB4D9",
  SQDG = "#9D4D6C",
  TG   = "#ED804A"
)

#---------------------------
# Valid classes (from paper script)
#---------------------------
valid_classes <- names(class_colors)
if (length(pca_class_choices) == 0) {
  pca_class_choices <- sort(valid_classes)
}

#---------------------------
# Function to extract lipid class from name (from paper script)
# e.g., "PC(34:1)" -> "PC", "TG(52:3)" -> "TG"
#---------------------------
get_lipid_class <- function(x) {
  cls <- stringr::str_extract(x, "^[A-Za-z0-9]+(?=\\()")
  cls <- ifelse(is.na(cls), NA_character_, cls)
  ifelse(cls %in% valid_classes, cls, "Other")
}

get_lipid_class_loose <- function(x) {
  cls <- stringr::str_extract(x, "^[A-Za-z0-9]+(?=\\()")
  cls[is.na(cls) | cls == ""] <- "Other"
  cls
}

normalize_lipid_key <- function(x) {
  tolower(gsub("[^a-z0-9]", "", as.character(x)))
}

build_lipid_class_maps <- function(class_df) {
  if (is.null(class_df) || nrow(class_df) == 0 || !all(c("Lipids", "Class") %in% names(class_df))) {
    return(list(exact = c(), normalized = c()))
  }

  exact_df <- class_df %>%
    dplyr::transmute(
      Lipids = as.character(Lipids),
      Class = as.character(Class)
    ) %>%
    dplyr::filter(!is.na(Lipids), Lipids != "", !is.na(Class), Class != "") %>%
    dplyr::distinct(Lipids, .keep_all = TRUE)

  norm_df <- exact_df %>%
    dplyr::mutate(norm_key = normalize_lipid_key(Lipids)) %>%
    dplyr::filter(norm_key != "") %>%
    dplyr::group_by(norm_key) %>%
    dplyr::summarise(Class = dplyr::first(Class), .groups = "drop")

  list(
    exact = stats::setNames(exact_df$Class, exact_df$Lipids),
    normalized = stats::setNames(norm_df$Class, norm_df$norm_key)
  )
}

lipid_class_maps <- build_lipid_class_maps(lipid_class_info)

#---------------------------
# Group colors (for ellipses - from paper script)
#---------------------------
group_colors <- c(
  "Glycolipids"    = "#E7298A",
  "Phospholipids"  = "#1B9E77",
  "Neutral"        = "#6A51A3",
  "Storage"        = "#ED804A"
)

#---------------------------
# Function to map class to group (from paper script)
#---------------------------
class_to_group <- function(cls) {
  dplyr::case_when(
    cls %in% c("PC", "PA", "PE", "PG", "PS", "Glycerophospholipid", "Sphingolipid", "Ether lipid") ~ "Phospholipids",
    cls %in% c("MGDG", "DGDG", "SQDG") ~ "Glycolipids",
    cls %in% c("DG", "MG", "Fatty acid", "Fatty acyl", "Prenol") ~ "Neutral",
    cls %in% c("TG", "Glycerolipid") ~ "Storage",
    TRUE ~ NA_character_
  )
}

infer_lipid_class <- function(x) {
  x_chr <- as.character(x)

  cls <- unname(lipid_class_maps$exact[x_chr])
  missing_idx <- which(is.na(cls) | cls == "")

  if (length(missing_idx) > 0) {
    norm_keys <- normalize_lipid_key(x_chr[missing_idx])
    cls_norm <- unname(lipid_class_maps$normalized[norm_keys])
    cls[missing_idx] <- cls_norm
  }

  missing_idx <- which(is.na(cls) | cls == "")
  if (length(missing_idx) > 0) {
    cls[missing_idx] <- get_lipid_class_loose(x_chr[missing_idx])
  }

  cls[is.na(cls) | cls == ""] <- "Other"
  cls
}

# Load gene annotation
gene_annotation <- read_gene_annotation(resolved_paths$gene_annotation)

# Load threshold-specific annotation RDS files (plog107/plog106/plog105)
annotation_threshold_choices <- c("7", "6", "5")
threshold_to_suffix <- c("7" = "107", "6" = "106", "5" = "105")

resolve_threshold_annotation_path <- function(dataset_name, trait_type, threshold) {
  suffix <- threshold_to_suffix[[as.character(threshold)]]
  if (is.null(suffix)) {
    return(NA_character_)
  }

  candidates <- if (trait_type == "individual") {
    c(
      file.path("data", "GWAS_RDS", "individual", paste0("all_annotations_", dataset_name, "_individual_plog", suffix, ".rds")),
      file.path("data", "GWAS_RDS", "individual", paste0("all_annotations_", dataset_name, "_individual.rds"))
    )
  } else {
    c(
      file.path("data", "GWAS_RDS", "sum_ratio", paste0("all_annotations_", dataset_name, "_sum_ratio_plog", suffix, ".rds")),
      file.path("data", "GWAS_RDS", "sum_ratio", paste0("all_annotations_", dataset_name, "_sum_ratio.rds"))
    )
  }

  resolve_existing_path(candidates)
}

build_annotation_store <- function(threshold) {
  lowinput_individual <- read_annotation_bundle(resolve_threshold_annotation_path("lowinput", "individual", threshold))
  lowinput_sum_ratio <- read_annotation_bundle(resolve_threshold_annotation_path("lowinput", "sum_ratio", threshold))
  control_individual <- read_annotation_bundle(resolve_threshold_annotation_path("control", "individual", threshold))
  control_sum_ratio <- read_annotation_bundle(resolve_threshold_annotation_path("control", "sum_ratio", threshold))

  lowinput_annotations <- merge_annotation_lists_with_source(
    individual = lowinput_individual,
    sum_ratio = lowinput_sum_ratio,
    other_label = "individual"
  )
  control_annotations <- merge_annotation_lists_with_source(
    individual = control_individual,
    sum_ratio = control_sum_ratio,
    other_label = "individual"
  )

  list(
    data = list(
      lowinput = lowinput_annotations$data,
      control = control_annotations$data
    ),
    source = list(
      lowinput = lowinput_annotations$source,
      control = control_annotations$source
    )
  )
}

annotation_store_by_threshold <- setNames(
  lapply(annotation_threshold_choices, build_annotation_store),
  annotation_threshold_choices
)

annotation_data_by_threshold <- lapply(annotation_store_by_threshold, function(x) x$data)
annotation_data_source_by_threshold <- lapply(annotation_store_by_threshold, function(x) x$source)

annotation_data_default <- annotation_data_by_threshold[["7"]]
if (is.null(annotation_data_default)) {
  annotation_data_default <- annotation_data_by_threshold[[annotation_threshold_choices[[1]]]]
}
if (is.null(annotation_data_default)) {
  annotation_data_default <- list(lowinput = list(), control = list())
}

annotation_dataset_choices <- names(annotation_data_default)[vapply(annotation_data_default, length, integer(1)) > 0]
if (length(annotation_dataset_choices) == 0) {
  annotation_dataset_choices <- names(annotation_data_default)
}

default_annotation_dataset <- if ("lowinput" %in% annotation_dataset_choices) {
  "lowinput"
} else {
  annotation_dataset_choices[[1]]
}

# Path to Manhattan figures
manhattan_dir <- file.path("www", "manhattans")
## If you're using prcomp, you can skip factoextra. If you want fancy visuals, use library(factoextra).
## If you're doing UMAP, library(umap) is also helpful.

## -- Load your data. Adjust paths if needed --
# e.g., lowinput <- vroom::vroom("../data/lowinput_exp.csv")
#       control  <- vroom::vroom("../data/control_exp.csv")

# For demonstration, we assume `control` and `lowinput` are in memory.

# Custom CSS for modern card styling
custom_css <- "
.setup-banner {
  max-width: 1400px;
  margin: 20px auto 0;
  padding: 18px 22px;
  border-radius: 14px;
  background: #fff3cd;
  border: 1px solid #f0d58c;
  color: #664d03;
  box-shadow: 0 6px 18px rgba(102, 77, 3, 0.08);
}

.setup-banner h4 {
  margin-top: 0;
  margin-bottom: 8px;
}

.setup-banner ul {
  margin-bottom: 12px;
}

/* Landing page styles */
.landing-container {
  padding: 40px 20px;
  max-width: 1400px;
  margin: 0 auto;
}

.landing-header {
  text-align: center;
  margin-bottom: 50px;
}

.landing-header h1 {
  font-size: 3rem;
  font-weight: 700;
  color: #2c3e50;
  margin-bottom: 10px;
}

.landing-header p {
  font-size: 1.2rem;
  color: #7f8c8d;
  max-width: 600px;
  margin: 0 auto;
}

.cards-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
  gap: 24px;
  padding: 20px 0;
}

.nav-card {
  background: white;
  border-radius: 16px;
  padding: 30px 24px;
  cursor: pointer;
  transition: all 0.3s ease;
  border: 1px solid #e9ecef;
  box-shadow: 0 2px 8px rgba(0,0,0,0.04);
  display: flex;
  flex-direction: column;
  align-items: center;
  text-align: center;
}

.nav-card:hover {
  transform: translateY(-6px);
  box-shadow: 0 12px 40px rgba(0,0,0,0.12);
  border-color: transparent;
}

.nav-card .card-icon {
  width: 70px;
  height: 70px;
  border-radius: 16px;
  display: flex;
  align-items: center;
  justify-content: center;
  margin-bottom: 20px;
  font-size: 28px;
  color: white;
}

.nav-card h3 {
  font-size: 1.25rem;
  font-weight: 600;
  color: #2c3e50;
  margin-bottom: 8px;
}

.nav-card p {
  font-size: 0.9rem;
  color: #95a5a6;
  line-height: 1.5;
  margin: 0;
}

/* Color themes for cards */
.card-intro .card-icon { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); }
.card-data .card-icon { background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%); }
.card-boxplot .card-icon { background: linear-gradient(135deg, #ee0979 0%, #ff6a00 100%); }
.card-correlation .card-icon { background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); }
.card-pca .card-icon { background: linear-gradient(135deg, #fa709a 0%, #fee140 100%); }
.card-tsne .card-icon { background: linear-gradient(135deg, #a8edea 0%, #fed6e3 100%); }
.card-tsne .card-icon { color: #666; }
.card-heatmap .card-icon { background: linear-gradient(135deg, #ff6b6b 0%, #feca57 100%); }
.card-gwas .card-icon { background: linear-gradient(135deg, #5f72bd 0%, #9b23ea 100%); }
.card-genehits .card-icon { background: linear-gradient(135deg, #36d1dc 0%, #5b86e5 100%); }
.card-genotype .card-icon { background: linear-gradient(135deg, #c471f5 0%, #fa71cd 100%); }
.card-venn .card-icon { background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); }
.card-volcano .card-icon { background: linear-gradient(135deg, #f5af19 0%, #f12711 100%); }
.card-network .card-icon { background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%); }

/* Legend for unique lipids */
.unique-lipid-legend {
  background: #f8f9fa;
  border-left: 4px solid #2d2d2d;
  padding: 10px 15px;
  margin: 10px 0;
  font-size: 13px;
  color: #555;
}
.unique-lipid-legend .legend-color {
  display: inline-block;
  width: 20px;
  height: 14px;
  background: #2d2d2d;
  vertical-align: middle;
  margin-right: 8px;
  border-radius: 2px;
}

/* Statistics badge */
.stat-badge {
  display: inline-block;
  padding: 4px 10px;
  border-radius: 12px;
  font-size: 12px;
  font-weight: 600;
  margin: 2px;
}
.stat-badge.significant { background: #d4edda; color: #155724; }
.stat-badge.not-significant { background: #f8d7da; color: #721c24; }

/* Page content styles */
.page-container {
  padding: 20px;
  max-width: 1600px;
  margin: 0 auto;
}

.page-header {
  display: flex;
  align-items: center;
  margin-bottom: 30px;
  padding-bottom: 20px;
  border-bottom: 1px solid #e9ecef;
}

.home-btn {
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  color: white !important;
  border: none;
  padding: 10px 20px;
  border-radius: 8px;
  font-weight: 500;
  cursor: pointer;
  display: flex;
  align-items: center;
  gap: 8px;
  transition: all 0.2s ease;
  text-decoration: none;
  margin-right: 20px;
}

.home-btn:hover {
  transform: scale(1.05);
  box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);
}

.page-title {
  font-size: 1.8rem;
  font-weight: 600;
  color: #2c3e50;
  margin: 0;
}

/* Modern card panels for content */
.content-card {
  background: white;
  border-radius: 12px;
  padding: 24px;
  margin-bottom: 20px;
  box-shadow: 0 2px 12px rgba(0,0,0,0.06);
  border: 1px solid #f0f0f0;
}

.content-card h4 {
  font-weight: 600;
  color: #2c3e50;
  margin-bottom: 16px;
  padding-bottom: 12px;
  border-bottom: 2px solid #f0f0f0;
}

/* Global body styling */
body {
  background-color: #f8f9fa;
  font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
}

/* Custom select and input styling */
.form-control, .selectize-input {
  border-radius: 8px !important;
  border: 1px solid #e0e0e0 !important;
}

.form-control:focus, .selectize-input.focus {
  border-color: #667eea !important;
  box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.1) !important;
}

/* Button styling */
.btn-primary {
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  border: none;
  border-radius: 8px;
}

.btn-primary:hover {
  background: linear-gradient(135deg, #5a6fd6 0%, #6a4190 100%);
}

/* Download button styling */
.download-section {
  display: flex;
  gap: 10px;
  margin-top: 15px;
  padding-top: 15px;
  border-top: 1px solid #e9ecef;
  flex-wrap: wrap;
}

.btn-download {
  background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
  color: white !important;
  border: none;
  padding: 8px 16px;
  border-radius: 8px;
  font-size: 0.85rem;
  font-weight: 500;
  cursor: pointer;
  display: inline-flex;
  align-items: center;
  gap: 6px;
  transition: all 0.2s ease;
  text-decoration: none;
}

.btn-download:hover {
  transform: scale(1.03);
  box-shadow: 0 4px 12px rgba(17, 153, 142, 0.3);
  color: white !important;
}

.btn-download-plot {
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
}

.btn-download-plot:hover {
  box-shadow: 0 4px 12px rgba(102, 126, 234, 0.3);
}

/* Top Navigation Bar */
.top-navbar {
  background: white;
  padding: 12px 20px;
  margin-bottom: 20px;
  border-radius: 12px;
  box-shadow: 0 2px 12px rgba(0,0,0,0.06);
  display: flex;
  align-items: center;
  gap: 8px;
  flex-wrap: wrap;
  border: 1px solid #f0f0f0;
}

.nav-item {
  padding: 8px 14px;
  border-radius: 8px;
  cursor: pointer;
  font-size: 0.85rem;
  font-weight: 500;
  color: #6c757d;
  transition: all 0.2s ease;
  display: flex;
  align-items: center;
  gap: 6px;
  text-decoration: none;
  white-space: nowrap;
}

.nav-item:hover {
  background: #f8f9fa;
  color: #2c3e50;
}

.nav-item.active {
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  color: white;
}

.nav-item.home-item {
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  color: white;
  margin-right: 8px;
}

.nav-item.home-item:hover {
  transform: scale(1.03);
  box-shadow: 0 4px 12px rgba(102, 126, 234, 0.3);
}

.nav-divider {
  width: 1px;
  height: 24px;
  background: #e9ecef;
  margin: 0 8px;
}

.page-title-section {
  margin-bottom: 24px;
}

.page-title-section h2 {
  font-size: 1.8rem;
  font-weight: 600;
  color: #2c3e50;
  margin: 0;
}

/* Contact page styling */
.contact-card {
  background: white;
  border-radius: 16px;
  padding: 40px;
  box-shadow: 0 4px 20px rgba(0,0,0,0.08);
  border: 1px solid #f0f0f0;
  max-width: 700px;
  margin: 0 auto;
}

.contact-header {
  text-align: center;
  margin-bottom: 30px;
}

.contact-header h3 {
  font-size: 1.8rem;
  font-weight: 600;
  color: #2c3e50;
  margin-bottom: 10px;
}

.contact-header p {
  color: #7f8c8d;
  font-size: 1rem;
}

.contact-info-item {
  display: flex;
  align-items: flex-start;
  gap: 16px;
  padding: 20px;
  background: #f8f9fa;
  border-radius: 12px;
  margin-bottom: 16px;
  transition: all 0.2s ease;
}

.contact-info-item:hover {
  background: #f0f4ff;
  transform: translateX(4px);
}

.contact-icon {
  width: 50px;
  height: 50px;
  border-radius: 12px;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 20px;
  color: white;
  flex-shrink: 0;
}

.contact-icon.location { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); }
.contact-icon.email { background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%); }
.contact-icon.phone { background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); }

.contact-details h4 {
  font-size: 1rem;
  font-weight: 600;
  color: #2c3e50;
  margin: 0 0 4px 0;
}

.contact-details p {
  color: #6c757d;
  margin: 0;
  font-size: 0.95rem;
  line-height: 1.5;
}

.contact-details a {
  color: #667eea;
  text-decoration: none;
  font-weight: 500;
}

.contact-details a:hover {
  text-decoration: underline;
}

.contact-footer {
  text-align: center;
  margin-top: 30px;
  padding-top: 20px;
  border-top: 1px solid #e9ecef;
}

.contact-footer p {
  color: #95a5a6;
  font-size: 0.9rem;
}

.card-contact .card-icon { background: linear-gradient(135deg, #ff6b6b 0%, #ee5a24 100%); }
"

# Helper function to create navigation cards
nav_card <- function(id, icon_name, title, description, card_class) {
  tags$div(
    class = paste("nav-card", card_class),
    onclick = sprintf("Shiny.setInputValue('nav_to', '%s', {priority: 'event'})", id),
    tags$div(class = "card-icon", icon(icon_name)),
    tags$h3(title),
    tags$p(description)
  )
}

# Helper function to create navigation item
nav_item <- function(id, icon_name, label, is_active = FALSE, is_home = FALSE) {
  class_name <- if (is_home) {
    "nav-item home-item"
  } else if (is_active) {
    "nav-item active"
  } else {
    "nav-item"
  }

  tags$a(
    class = class_name,
    onclick = sprintf("Shiny.setInputValue('nav_to', '%s', {priority: 'event'})", id),
    icon(icon_name),
    label
  )
}

# Helper function to create page header with full navigation bar
page_header <- function(title, current_tab = "") {
  tags$div(
    # Navigation bar
    tags$div(
      class = "top-navbar",
      nav_item("home", "home", "Home", is_home = TRUE),
      tags$div(class = "nav-divider"),
      nav_item("intro", "info-circle", "Intro", current_tab == "intro"),
      nav_item("data", "table", "Data", current_tab == "data"),
      nav_item("boxplot", "chart-bar", "Boxplots", current_tab == "boxplot"),
      nav_item("correlation", "project-diagram", "Correlation", current_tab == "correlation"),
      nav_item("pca", "compress", "PCA", current_tab == "pca"),
      nav_item("tsne", "braille", "t-SNE", current_tab == "tsne"),
      nav_item("heatmap", "th", "Heatmap", current_tab == "heatmap"),
      nav_item("gwas", "dna", "GWAS", current_tab == "gwas"),
      nav_item("genehits", "bullseye", "Gene Hits", current_tab == "genehits"),
      nav_item("genotype", "filter", "Genotype", current_tab == "genotype"),
      nav_item("venn", "object-group", "Venn", current_tab == "venn"),
      nav_item("volcano", "mountain", "Volcano", current_tab == "volcano"),
      nav_item("network", "share-alt", "Network", current_tab == "network"),
      nav_item("contact", "envelope", "Contact", current_tab == "contact")
    ),
    # Page title
    tags$div(
      class = "page-title-section",
      tags$h2(title)
    )
  )
}

# Helper function for content cards
content_card <- function(..., title = NULL) {
  tags$div(
    class = "content-card",
    if (!is.null(title)) tags$h4(title),
    ...
  )
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML(custom_css)),
    tags$link(
      href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap",
      rel = "stylesheet"
    )
  ),
  uiOutput("missing_data_banner"),

  # Hidden navigation tabs
  tabsetPanel(
    id = "main_tabs",
    type = "hidden",

    # ============================================
    # LANDING PAGE
    # ============================================
    tabPanel(
      "home",
      tags$div(
        class = "landing-container",
        tags$div(
          class = "landing-header",
          tags$h1(app_title),
          tags$p("Explore SAP sorghum lipid profiles, GWAS hits, and comparative plots across field conditions")
        ),
        tags$div(
          class = "cards-grid",
          nav_card("intro", "info-circle", "Introduction",
                   "Learn about the experiment and dataset", "card-intro"),
          nav_card("data", "table", "Data Preview",
                   "Browse and filter lipid abundance data", "card-data"),
          nav_card("boxplot", "chart-bar", "Boxplots",
                   "Visualize lipid distributions", "card-boxplot"),
          nav_card("correlation", "project-diagram", "Correlation",
                   "Analyze lipid correlations", "card-correlation"),
          nav_card("pca", "compress", "PCA",
                   "Principal component analysis", "card-pca"),
          nav_card("tsne", "braille", "t-SNE / UMAP",
                   "Dimensionality reduction plots", "card-tsne"),
          nav_card("heatmap", "th", "Heatmap",
                   "Visualize expression patterns", "card-heatmap"),
          nav_card("gwas", "dna", "GWAS",
                   "Genome-wide association results", "card-gwas"),
          nav_card("genehits", "bullseye", "Gene Hits",
                   "Top genes by association count", "card-genehits"),
          nav_card("genotype", "filter", "Genotype",
                   "Filter genotypes by lipid levels", "card-genotype"),
          nav_card("venn", "object-group", "Venn Diagram",
                   "Compare lipids across conditions", "card-venn"),
          nav_card("volcano", "mountain", "Volcano Plot",
                   "Differential lipid analysis", "card-volcano"),
          nav_card("network", "share-alt", "Network",
                   "Lipid correlation networks", "card-network"),
          nav_card("contact", "envelope", "Contact Us",
                   "Get in touch with our team", "card-contact")
        )
      )
    ),
      
      
    # ============================================
    # INTRODUCTION PAGE
    # ============================================
    tabPanel(
      "intro",
      tags$div(
        class = "page-container",
        page_header("Introduction", "intro"),
        fluidRow(
          column(12,
            content_card(
              title = paste("Welcome to the", app_title),
              p("This Shiny app presents lipidomics data generated from an experiment using the Sorghum Association Panel (SAP), a genetically diverse set of sorghum lines."),
              p("Our goal is to provide a reference lipid profile under two distinct field conditions, enabling researchers to explore, compare, and utilize lipid data from sorghum grown under contrasting environments."),
              br(),
              h4("Experimental Setup"),
              p("We planted SAP lines in two field conditions over two growing seasons (2019 and 2022):"),
              tags$ul(
                tags$li(strong("Control Condition:"), " Standard agricultural input with adequate nitrogen (N), phosphorus (P), and normal planting time."),
                tags$li(strong("Low Input Condition:"), " Characterized by low nitrogen, low phosphorus, and early planting to mimic cooler, more stressful conditions.")
              ),
              br(),
              p("After the plants matured, we performed lipidomics analysis, profiling more than 250 lipid compounds from each condition. These include:"),
              tags$ul(
                tags$li("Condition-specific lipids unique to either control or low-input fields"),
                tags$li("Shared lipids found in both conditions but varying in abundance")
              ),
              br(),
              h4("What You Can Do with This App"),
              tags$ul(
                tags$li("Explore lipid profiles across different genotypes and conditions"),
                tags$li("Visualize condition-specific lipid signatures"),
                tags$li("Download curated datasets"),
                tags$li("Use the data as a reference for your own sorghum lipidomics experiments")
              ),
              br(),
              p("Our long-term goal is to establish a publicly accessible sorghum lipid database to support the research community in advancing sorghum studies, and metabolic studies.")
            )
          )
        )
      )
    ),
      
      
    # ============================================
    # DATA PREVIEW PAGE
    # ============================================
    tabPanel(
      "data",
      tags$div(
        class = "page-container",
        page_header("Data Preview", "data"),
        tabsetPanel(
          type = "pills",
          tabPanel(
            "Data Table",
            br(),
            fluidRow(
              column(4,
                content_card(
                  title = "Choose Dataset",
                  selectInput("data_choice", "Dataset:", choices = c("control", "lowinput")),
                  selectInput("data_type_data", "Data Type:",
                              choices = data_type_choices,
                              selected = "raw"),
                  uiOutput("class_filter_ui"),
                  uiOutput("subclass_filter_data_ui"),
                  uiOutput("lipid_filter_data_ui")
                )
              ),
              column(8,
                content_card(
                  title = "Table Preview",
                  DTOutput("data_table"),
                  tags$div(
                    class = "download-section",
                    downloadButton("download_data_table", "Download CSV", class = "btn-download")
                  )
                )
              )
            )
          ),
          tabPanel(
            "Data Visualize",
            br(),
            fluidRow(
              column(4,
                content_card(
                  title = "Visualization Controls",
                  selectInput("data_viz_domain", "Data Domain:",
                              choices = c("Lipid Abundance" = "abundance", "GWAS Traits" = "gwas"),
                              selected = "abundance"),
                  selectInput("data_viz_mode", "View Mode:",
                              choices = c("Single Dataset" = "single", "Compare Control vs Low Input" = "compare"),
                              selected = "single"),
                  conditionalPanel(
                    condition = "input.data_viz_mode == 'single'",
                    selectInput("data_viz_dataset", "Dataset:", choices = c("control", "lowinput"), selected = "control")
                  ),
                  conditionalPanel(
                    condition = "input.data_viz_domain == 'abundance'",
                    selectInput("data_viz_data_type", "Data Type:",
                                choices = data_type_choices,
                                selected = "tic")
                  ),
                  conditionalPanel(
                    condition = "input.data_viz_domain == 'gwas'",
                    selectInput("data_viz_gwas_threshold", "GWAS -log10(p) threshold:",
                                choices = c("7", "6", "5"),
                                selected = "7"),
                    selectInput("data_viz_trait_source", "Trait Source:",
                                choices = c("All" = "all", "Individual" = "individual", "Sum/Ratio" = "sum_ratio"),
                                selected = "all")
                  ),
                  selectInput("data_viz_level", "Visualize By:",
                              choices = c("Class", "Subclass", "Lipid"),
                              selected = "Class"),
                  uiOutput("data_viz_class_filter_ui"),
                  uiOutput("data_viz_subclass_filter_ui"),
                  uiOutput("data_viz_lipid_filter_ui"),
                  sliderInput("data_viz_top_n", "Top groups for plot:", min = 5, max = 30, value = 12, step = 1),
                  helpText("Single mode draws one stacked bar. Compare mode draws Control and Low Input bars.")
                )
              ),
              column(8,
                content_card(
                  title = "Data Visualization",
                  withSpinner(plotOutput("data_visualize_plot", height = "520px"), type = 6, color = "#3498db"),
                  uiOutput("data_visualize_summary")
                )
              )
            )
          )
        )
      )
    ),

    # ============================================
    # BOXPLOTS PAGE
    # ============================================
    tabPanel(
      "boxplot",
      tags$div(
        class = "page-container",
        page_header("Boxplots", "boxplot"),
        tabsetPanel(
          type = "pills",
          tabPanel("Absolute",
            br(),
            fluidRow(
              column(4,
                content_card(
                  title = "Plot Controls",
                  selectInput("plot_dataset", "Choose Dataset:", choices = c("control", "lowinput")),
                  selectInput("data_type_boxplot", "Data Type:",
                              choices = data_type_choices,
                              selected = "raw"
                  ),
                  uiOutput("class_filter_plot_ui"),
                  uiOutput("subclass_filter_ui"),
                  uiOutput("subsubclass_filter_ui"),
                  uiOutput("compound_ui")
                )
              ),
              column(8,
                content_card(
                  title = "Distribution Plots",
                  plotOutput("boxplot_single"),
                  plotOutput("boxplot_both"),
                  tags$div(
                    class = "download-section",
                    downloadButton("download_boxplot_single", "Download Density Plot", class = "btn-download btn-download-plot"),
                    downloadButton("download_boxplot_both", "Download Comparison Plot", class = "btn-download btn-download-plot")
                  )
                )
              )
            )
          ),
          tabPanel("Relative",
            br(),
            fluidRow(
              column(4,
                content_card(
                  title = "Ratio Controls",
                  selectInput(
                    inputId = "plot_dataset_ratio",
                    label = "Choose Dataset:",
                    choices = c("control", "lowinput")
                  ),
                  uiOutput("ratio_ui_1"),
                  uiOutput("ratio_ui_2")
                )
              ),
              column(8,
                content_card(
                  title = "Ratio Boxplots",
                  plotOutput("boxplot_ratio"),
                  plotOutput("boxplot_ratio_combined"),
                  tags$div(
                    class = "download-section",
                    downloadButton("download_ratio_single", "Download Ratio Plot", class = "btn-download btn-download-plot"),
                    downloadButton("download_ratio_combined", "Download Combined Plot", class = "btn-download btn-download-plot")
                  )
                )
              )
            )
          )
        )
      )
    ),


    # ============================================
    # PCA PAGE
    # ============================================
    tabPanel(
      "pca",
      tags$div(
        class = "page-container",
        page_header("Principal Component Analysis", "pca"),
        fluidRow(
          column(4,
            content_card(
              title = "PCA Controls",
              selectInput(
                inputId = "pca_dataset",
                label = "Select Dataset:",
                choices = c("control", "lowinput")
              ),
              selectInput(
                inputId = "data_type_pca",
                label = "Data Type:",
                choices = data_type_choices,
                selected = "raw"
              ),
              checkboxInput("pca_center", "Center", TRUE),
              checkboxInput("pca_scale", "Scale", TRUE),
              selectInput(
                inputId = "pca_mode",
                label = "PCA Mode:",
                choices = c("Samples", "Lipids"),
                selected = "Samples"
              ),
              conditionalPanel(
                condition = "input.pca_mode == 'Lipids'",
                selectInput(
                  inputId = "pca_class_filter",
                  label = "Filter by Lipid Class:",
                  choices = c("All", pca_class_choices),
                  selected = "All"
                ),
                helpText("Points colored by class, ellipses by group (Phospholipids, Glycolipids, Neutral, Storage)")
              ),
              helpText("Raw and %TIC use standard PCA controls."),
              helpText("CLR applies closure + pseudocount + CLR + z-scoring before PCA."),
              verbatimTextOutput("pca_selected")
            )
          ),
          column(8,
            content_card(
              title = "PCA Plot",
              withSpinner(plotlyOutput("pca_plot"), type = 6, color = "#3498db"),
              uiOutput("pca_variance_info"),
              tags$div(
                class = "download-section",
                downloadButton("download_pca_plot", "Download PCA Plot (PNG)", class = "btn-download btn-download-plot")
              )
            )
          )
        )
      )
    ),



    # ============================================
    # t-SNE / UMAP PAGE
    # ============================================
    tabPanel(
      "tsne",
      tags$div(
        class = "page-container",
        page_header("t-SNE / UMAP", "tsne"),
        fluidRow(
          column(4,
            content_card(
              title = "Dimensionality Reduction Controls",
              selectInput(
                inputId = "dr_dataset",
                label = "Select Dataset:",
                choices = c("control", "lowinput")
              ),
              selectInput(
                inputId = "data_type_dr",
                label = "Data Type:",
                choices = data_type_choices,
                selected = "raw"
              ),
              selectInput(
                inputId = "dr_method",
                label = "Method:",
                choices = c("t-SNE", "UMAP")
              ),
              selectInput(
                inputId = "dr_mode",
                label = "Dimension Mode:",
                choices = c("Samples", "Lipids"),
                selected = "Samples"
              ),
              numericInput("perplexity", "t-SNE perplexity:", 10, min = 1, max = 50),
              helpText("For t-SNE, perplexity is used. For UMAP, it's ignored."),
              conditionalPanel(
                condition = "input.dr_mode == 'Lipids'",
                helpText("Points colored by class, ellipses by group (Phospholipids, Glycolipids, Neutral, Storage)")
              ),
              verbatimTextOutput("dr_selected")
            )
          ),
          column(8,
            content_card(
              title = "Visualization",
              withSpinner(plotlyOutput("dr_plot"), type = 6, color = "#3498db"),
              tags$div(
                class = "download-section",
                downloadButton("download_dr_plot", "Download Plot (PNG)", class = "btn-download btn-download-plot")
              )
            )
          )
        )
      )
    ),



    # ============================================
    # HEATMAP PAGE
    # ============================================
    tabPanel(
      "heatmap",
      tags$div(
        class = "page-container",
        page_header("Heatmap Viewer", "heatmap"),
        tabsetPanel(
          type = "pills",

          tabPanel("Individual Samples",
            br(),
            fluidRow(
              column(4,
                content_card(
                  title = "Heatmap Controls",
                  selectInput("heatmap_dataset", "Select Dataset:", choices = c("control", "lowinput", "Both")),
                  selectInput("data_type_heatmap", "Data Type:",
                              choices = data_type_choices,
                              selected = "raw"),
                  selectInput("heatmap_mode", "With Respect To:", choices = c("Samples", "Lipids"), selected = "Samples"),
                  selectInput("heatmap_lipid_selection", "Lipid Selection:",
                              choices = c(
                                "Top Changed Lipids" = "top_changed",
                                "By Lipid Class" = "by_class",
                                "Individual Lipids" = "individual"
                              ),
                              selected = "top_changed"),
                  conditionalPanel(
                    condition = "input.heatmap_lipid_selection == 'top_changed'",
                    sliderInput("top_n_heatmap", "Top Changed Lipids:", min = 5, max = 100, value = 20),
                    helpText("Shows top changed lipids based on mean absolute deviation from median.")
                  ),
                  conditionalPanel(
                    condition = "input.heatmap_lipid_selection == 'by_class'",
                    selectInput("heatmap_class_filter", "Select Lipid Class:",
                                choices = sort(valid_classes),
                                selected = "PC"),
                    helpText("Shows all lipids from the selected class.")
                  ),
                  conditionalPanel(
                    condition = "input.heatmap_lipid_selection == 'individual'",
                    uiOutput("heatmap_lipid_selector"),
                    helpText("Select specific lipids to display in the heatmap.")
                  ),
                  selectInput("heatmap_color_palette", "Color Palette:",
                              choices = c("Red-Yellow-Blue" = "RdYlBu",
                                          "Viridis (Colorblind-friendly)" = "viridis",
                                          "Plasma" = "plasma"),
                              selected = "RdYlBu")
                )
              ),
              column(8,
                content_card(
                  title = "Heatmap (Samples x Lipids)",
                  conditionalPanel(
                    condition = "input.heatmap_dataset == 'Both'",
                    tags$div(class = "unique-lipid-legend",
                      tags$span(class = "legend-color"),
                      "Dark gray cells indicate lipids not present in that condition"
                    )
                  ),
                  withSpinner(plotlyOutput("heatmap_plot", height = "700px"), type = 6, color = "#3498db"),
                  DT::dataTableOutput("lipid_name_lookup"),
                  tags$div(
                    class = "download-section",
                    downloadButton("download_lipid_lookup", "Download Lipid Key (CSV)", class = "btn-download")
                  )
                )
              )
            )
          ),

          tabPanel("Condition Comparison",
            br(),
            fluidRow(
              column(4,
                content_card(
                  title = "Comparison Controls",
                  selectInput("data_type_heatmap_diff", "Data Type:",
                              choices = data_type_choices,
                              selected = "raw"),
                  selectInput("heatmap_diff_lipid_selection", "Lipid Selection:",
                              choices = c(
                                "Top Changed Lipids" = "top_changed",
                                "By Lipid Class" = "by_class",
                                "Individual Lipids" = "individual"
                              ),
                              selected = "top_changed"),
                  conditionalPanel(
                    condition = "input.heatmap_diff_lipid_selection == 'top_changed'",
                    sliderInput("top_n_diff", "Top Differential Lipids:", min = 5, max = 100, value = 20),
                    helpText("Lipids with largest Z-score difference between control and lowinput.")
                  ),
                  conditionalPanel(
                    condition = "input.heatmap_diff_lipid_selection == 'by_class'",
                    selectInput("heatmap_diff_class_filter", "Select Lipid Class:",
                                choices = sort(valid_classes),
                                selected = "PC"),
                    helpText("Shows all lipids from the selected class.")
                  ),
                  conditionalPanel(
                    condition = "input.heatmap_diff_lipid_selection == 'individual'",
                    uiOutput("heatmap_diff_lipid_selector"),
                    helpText("Select specific lipids to display in the heatmap.")
                  ),
                  selectInput("heatmap_diff_color_palette", "Color Palette:",
                              choices = c("Red-Yellow-Blue" = "RdYlBu",
                                          "Viridis (Colorblind-friendly)" = "viridis",
                                          "Plasma" = "plasma"),
                              selected = "RdYlBu")
                )
              ),
              column(8,
                content_card(
                  title = "Comparison Heatmap",
                  tags$div(class = "unique-lipid-legend",
                    tags$span(class = "legend-color"),
                    "Dark gray cells indicate lipids unique to one condition (not present in the other)"
                  ),
                  withSpinner(plotlyOutput("heatmap_diff", height = "700px"), type = 6, color = "#3498db"),
                  DT::dataTableOutput("lipid_name_lookup_diff"),
                  tags$div(
                    class = "download-section",
                    downloadButton("download_lipid_lookup_diff", "Download Lipid Key (CSV)", class = "btn-download")
                  )
                )
              )
            )
          ),

          tabPanel("Summed Classes",
            br(),
            fluidRow(
              column(4,
                content_card(
                  title = "Summed Class Controls",
                  selectInput("heatmap_summed_dataset", "Select Dataset:", choices = c("control", "lowinput")),
                  selectInput("data_type_heatmap_summed", "Data Type:",
                              choices = data_type_choices,
                              selected = "raw"),
                  selectInput("heatmap_summed_mode", "With Respect To:", choices = c("Samples", "Lipids"), selected = "Samples"),
                  sliderInput("top_n_summed", "Top Changed Classes:", min = 5, max = 30, value = 10),
                  helpText("Each column represents total lipid abundance per class.")
                )
              ),
              column(8,
                content_card(
                  title = "Summed Lipid Class Heatmap",
                  withSpinner(plotlyOutput("heatmap_summed", height = "700px"), type = 6, color = "#3498db")
                )
              )
            )
          )
        )
      )
    ),

    # ============================================
    # GWAS PAGE
    # ============================================
    tabPanel(
      "gwas",
      tags$div(
        class = "page-container",
        page_header("GWAS Results", "gwas"),
        fluidRow(
          column(4,
            content_card(
              title = "GWAS Controls",
              selectInput("gwas_dataset", "Select Dataset:", choices = annotation_dataset_choices, selected = default_annotation_dataset),
              selectInput("gwas_threshold", "Select -log10(p) threshold:", choices = c("7", "6", "5"), selected = "7"),
              selectInput("gwas_trait_source", "Trait Source:", choices = c("All" = "all", "Individual" = "individual", "Sum/Ratio" = "sum_ratio"), selected = "all"),
              uiOutput("gwas_class_filter"),
              uiOutput("gwas_subclass_filter"),
              uiOutput("gwas_lipid_selector")
            )
          ),
          column(8,
            content_card(
              title = "Manhattan Plot",
              uiOutput("manhattan_image")
            )
          )
        ),
        fluidRow(
          column(12,
            content_card(
              title = "GWAS Annotation Table",
              uiOutput("gwas_trait_class_info"),
              DTOutput("gwas_table"),
              tags$div(
                class = "download-section",
                downloadButton("download_gwas_table", "Download GWAS Results (CSV)", class = "btn-download")
              )
            )
          )
        )
      )
    ),

    # ============================================
    # GENE HITS PAGE
    # ============================================
    tabPanel(
      "genehits",
      tags$div(
        class = "page-container",
        page_header("Gene Hits", "genehits"),
        fluidRow(
          column(4,
            content_card(
              title = "Gene Hit Filters",
              selectInput("hit_dataset", "Select Dataset:", choices = annotation_dataset_choices, selected = default_annotation_dataset),
              selectInput("hit_trait_source", "Trait Source:", choices = c("All" = "all", "Individual" = "individual", "Sum/Ratio" = "sum_ratio"), selected = "all"),
              selectInput("hit_threshold", "Select -log10(p) threshold:", choices = c("7", "6", "5"), selected = "7"),
              uiOutput("hit_class_filter"),
              uiOutput("hit_subclass_filter"),
              uiOutput("hit_lipid_filter")
            )
          ),
          column(8,
            content_card(
              title = "Gene and Phenotype Summaries",
              tabsetPanel(
                type = "pills",
                tabPanel(
                  "Top Genes",
                  br(),
                  DTOutput("gene_hit_table"),
                  tags$div(
                    class = "download-section",
                    downloadButton("download_gene_hits", "Download Gene Hits (CSV)", class = "btn-download")
                  )
                ),
                tabPanel(
                  "Phenotypes",
                  br(),
                  DTOutput("phenotype_hit_table"),
                  tags$div(
                    class = "download-section",
                    downloadButton("download_phenotype_hits", "Download Phenotypes (CSV)", class = "btn-download")
                  )
                )
              )
            )
          )
        )
      )
    ),

    # ============================================
    # CORRELATION PAGE
    # ============================================
    tabPanel(
      "correlation",
      tags$div(
        class = "page-container",
        page_header("Correlation Analysis", "correlation"),
        tabsetPanel(
          type = "pills",

          tabPanel("Individual Lipids",
            br(),
            fluidRow(
              column(4,
                content_card(
                  title = "Correlation Controls",
                  selectInput("cor_dataset", "Select Dataset:", choices = names(dataset_list), selected = "lowinput"),
                  selectInput("data_type_cor", "Data Type:",
                              choices = data_type_choices,
                              selected = "raw"),
                  selectInput("cor_class_filter", "Filter by Lipid Class:",
                              choices = c("All Classes" = "all", sort(valid_classes)),
                              selected = "all"),
                  uiOutput("lipid1_selector"),
                  uiOutput("lipid2_selector"),
                  helpText("Filter lipids by class, or select 'All Classes' to see all lipids.")
                )
              ),
              column(8,
                content_card(
                  title = "Lipid Correlation Plot",
                  plotOutput("correlation_plot"),
                  tags$div(
                    class = "download-section",
                    downloadButton("download_correlation_plot", "Download Correlation Plot", class = "btn-download btn-download-plot")
                  )
                )
              )
            )
          ),

          tabPanel("Summed Classes",
            br(),
            fluidRow(
              column(12,
                content_card(
                  title = "Class-wise Correlation Heatmap",
                  selectInput("class_cor_dataset", "Select Dataset:", choices = names(dataset_list)),
                  selectInput("data_type_class_cor", "Data Type:",
                              choices = data_type_choices,
                              selected = "raw"),
                  withSpinner(plotlyOutput("class_correlation_heatmap"), type = 6, color = "#3498db")
                )
              )
            )
          ),

          tabPanel("Correlation Heatmap",
            br(),
            fluidRow(
              column(4,
                content_card(
                  title = "Top Correlated Lipids",
                  selectInput("corheat_dataset", "Select Dataset:", choices = names(dataset_list)),
                  selectInput("data_type_corheat", "Data Type:",
                              choices = data_type_choices,
                              selected = "raw"),
                  sliderInput("top_mad_lipids", "Top Lipids (based on MAD):", min = 5, max = 100, value = 30)
                )
              ),
              column(8,
                content_card(
                  title = "Correlation Heatmap",
                  withSpinner(plotlyOutput("individual_correlation_heatmap"), type = 6, color = "#3498db")
                )
              )
            )
          )
        )
      )
    ),

    # ============================================
    # GENOTYPE PAGE
    # ============================================
    tabPanel(
      "genotype",
      tags$div(
        class = "page-container",
        page_header("Genotype Analysis", "genotype"),
        tabsetPanel(
          type = "pills",

          tabPanel("Genotype Selection",
            br(),
            fluidRow(
              column(4,
                content_card(
                  title = "Selection Controls",
                  uiOutput("geno_dataset_ui"),
                  selectInput("data_type_geno", "Data Type:",
                              choices = data_type_choices,
                              selected = "raw"),
                  uiOutput("geno_lipid_ui"),
                  selectInput("geno_filter_method", "Filter Type:",
                              choices = c("Custom Range", "Lowest Quartile", "Highest Quartile")),
                  conditionalPanel(
                    condition = "input.geno_filter_method == 'Custom Range'",
                    uiOutput("geno_slider_ui")
                  ),
                  actionButton("filter_genotypes", "Submit", class = "btn-primary", icon = icon("check"))
                )
              ),
              column(8,
                content_card(
                  title = "Matching Genotypes",
                  verbatimTextOutput("genotype_output")
                )
              )
            )
          ),

          tabPanel("Genotype Viewer",
            br(),
            fluidRow(
              column(4,
                content_card(
                  title = "View Genotype Lipids",
                  uiOutput("geno_dataset_ui2"),
                  selectInput("data_type_geno_view", "Data Type:",
                              choices = data_type_choices,
                              selected = "raw"),
                  uiOutput("geno_viewer_ui"),
                  uiOutput("geno_view_lipid_ui"),
                  selectInput("geno_ref_lines", "Show Reference Lines:",
                              choices = c(
                                "None" = "none",
                                "Mean" = "mean",
                                "Median" = "median",
                                "1st Quartile (Q1)" = "q1",
                                "3rd Quartile (Q3)" = "q3",
                                "All Quartiles" = "all"
                              ),
                              selected = "none"),
                  helpText("Reference lines show statistics across all genotypes for the selected lipid.")
                )
              ),
              column(8,
                content_card(
                  title = "Lipid Level Across Conditions",
                  plotOutput("geno_view_plot"),
                  tags$div(
                    class = "download-section",
                    downloadButton("download_geno_plot", "Download Plot", class = "btn-download btn-download-plot")
                  )
                )
              )
            )
          )
        )
      )
    ),

    # ============================================
    # VENN DIAGRAM PAGE
    # ============================================
    tabPanel(
      "venn",
      tags$div(
        class = "page-container",
        page_header("Venn Diagram", "venn"),
        fluidRow(
          column(4,
            content_card(
              title = "Venn Controls",
              uiOutput("venn_lipid_type_ui"),
              conditionalPanel(
                condition = "input.venn_type == 'Class'",
                uiOutput("venn_class_selector")
              ),
              helpText("Compare lipid compounds or lipid classes between Control and Low Input.")
            )
          ),
          column(8,
            content_card(
              title = "Venn Diagram",
              plotOutput("venn_plot"),
              tags$div(
                class = "download-section",
                downloadButton("download_venn_plot", "Download Venn Diagram", class = "btn-download btn-download-plot")
              )
            )
          )
        ),
        fluidRow(
          column(12,
            content_card(
              title = "Overlapping and Unique Sets",
              DTOutput("venn_table"),
              tags$div(
                class = "download-section",
                downloadButton("download_venn_table", "Download Venn Sets (CSV)", class = "btn-download"),
                downloadButton("download_combined_data", "Download Combined Lipid List (CSV)", class = "btn-download")
              )
            )
          )
        )
      )
    ),

    # ============================================
    # VOLCANO PLOT PAGE
    # ============================================
    tabPanel(
      "volcano",
      tags$div(
        class = "page-container",
        page_header("Volcano Plot", "volcano"),
        fluidRow(
          column(4,
            content_card(
              title = "Volcano Plot Controls",
              selectInput("volcano_data_type", "Data Type:",
                          choices = data_type_choices,
                          selected = "raw"),
              sliderInput("volcano_fc_threshold", "Log2 Fold Change Threshold:",
                          min = 0.5, max = 3, value = 1, step = 0.25),
              sliderInput("volcano_pval_threshold", "-log10(p-value) Threshold:",
                          min = 1, max = 5, value = 1.3, step = 0.1),
              checkboxInput("volcano_show_labels", "Show Lipid Labels", FALSE),
              sliderInput("volcano_top_labels", "Number of Labels to Show:",
                          min = 5, max = 50, value = 10),
              selectInput("volcano_color_mode", "Color Mode:",
                          choices = c("Standard" = "standard", "Colorblind-friendly" = "viridis"),
                          selected = "standard"),
              helpText("Compares Low Input vs Control. Positive FC = higher in Low Input.")
            )
          ),
          column(8,
            content_card(
              title = "Differential Lipid Analysis",
              withSpinner(plotlyOutput("volcano_plot", height = "600px"), type = 6, color = "#3498db"),
              uiOutput("volcano_summary"),
              tags$div(
                class = "download-section",
                downloadButton("download_volcano_plot", "Download Plot (PNG)", class = "btn-download btn-download-plot"),
                downloadButton("download_volcano_data", "Download Results (CSV)", class = "btn-download")
              )
            )
          )
        ),
        fluidRow(
          column(12,
            content_card(
              title = "Differential Lipids Table",
              DTOutput("volcano_table")
            )
          )
        ),
        fluidRow(
          column(12,
            content_card(
              title = "Lipid Class Enrichment Analysis",
              withSpinner(plotOutput("enrichment_plot", height = "400px"), type = 6, color = "#3498db"),
              helpText("Shows whether certain lipid classes are over/under-represented among significantly changed lipids.")
            )
          )
        )
      )
    ),

    # ============================================
    # NETWORK PAGE
    # ============================================
    tabPanel(
      "network",
      tags$div(
        class = "page-container",
        page_header("Lipid Correlation Network", "network"),
        fluidRow(
          column(4,
            content_card(
              title = "Network Controls",
              selectInput("network_dataset", "Select Dataset:",
                          choices = c("control", "lowinput")),
              selectInput("network_data_type", "Data Type:",
                          choices = data_type_choices,
                          selected = "raw"),
              selectInput("network_method", "Correlation Method:",
                          choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
                          selected = "pearson"),
              sliderInput("network_cor_threshold", "Correlation Threshold:",
                          min = 0.5, max = 0.95, value = 0.7, step = 0.05),
              selectInput("network_lipid_selection", "Lipid Selection:",
                          choices = c("Top Variable Lipids" = "top_var",
                                      "By Lipid Class" = "by_class"),
                          selected = "top_var"),
              conditionalPanel(
                condition = "input.network_lipid_selection == 'top_var'",
                sliderInput("network_top_n", "Number of Lipids:", min = 20, max = 100, value = 50)
              ),
              conditionalPanel(
                condition = "input.network_lipid_selection == 'by_class'",
                selectInput("network_class_filter", "Select Lipid Class:",
                            choices = sort(valid_classes), selected = "PC")
              ),
              checkboxInput("network_show_negative", "Show Negative Correlations", FALSE),
              helpText("Nodes are lipids, edges connect correlated lipids. Node size = number of connections.")
            )
          ),
          column(8,
            content_card(
              title = "Correlation Network",
              withSpinner(visNetworkOutput("network_plot", height = "600px"), type = 6, color = "#3498db"),
              tags$div(
                class = "download-section",
                downloadButton("download_network_edges", "Download Edge List (CSV)", class = "btn-download")
              )
            )
          )
        ),
        fluidRow(
          column(6,
            content_card(
              title = "Network Statistics",
              verbatimTextOutput("network_stats")
            )
          ),
          column(6,
            content_card(
              title = "Hub Lipids (Most Connected)",
              DTOutput("network_hubs")
            )
          )
        )
      )
    ),

    # ============================================
    # CONTACT US PAGE
    # ============================================
    tabPanel(
      "contact",
      tags$div(
        class = "page-container",
        page_header("Contact Us", "contact"),
        fluidRow(
          column(12,
            tags$div(
              class = "contact-card",
              tags$div(
                class = "contact-header",
                tags$h3("Get in Touch"),
                tags$p("Have questions about the data, methodology, or potential collaborations? We'd love to hear from you!")
              ),

              # Location
              tags$div(
                class = "contact-info-item",
                tags$div(class = "contact-icon location", icon("location-dot")),
                tags$div(
                  class = "contact-details",
                  tags$h4("Location"),
                  tags$p("North Carolina State University"),
                  tags$p("Department of Structural and Molecular Biochemistry"),
                  tags$p("Polk Hall, Room 134, Room 145 (Lab), Room 148 (Office)")
                )
              ),

              # Email
              tags$div(
                class = "contact-info-item",
                tags$div(class = "contact-icon email", icon("envelope")),
                tags$div(
                  class = "contact-details",
                  tags$h4("Email"),
                  tags$p(tags$a(href = "mailto:rrellan@ncsu.edu", "rrellan@ncsu.edu"))
                )
              ),

              # Phone
              tags$div(
                class = "contact-info-item",
                tags$div(class = "contact-icon phone", icon("phone")),
                tags$div(
                  class = "contact-details",
                  tags$h4("Phone"),
                  tags$p(tags$a(href = "tel:+19195154738", "(919) 515-4738"))
                )
              ),

              # Footer
              tags$div(
                class = "contact-footer",
                tags$p("We typically respond within 1-2 business days."),
                tags$p("For data access requests or collaboration inquiries, please include details about your research in your message.")
              )
            )
          )
        )
      )
    )

  )
)

server <- function(input, output, session) {

  output$missing_data_banner <- renderUI({
    build_setup_banner()
  })

  ##----------------------------------------------------------------------
  ##  Navigation handler for landing page cards
  ##----------------------------------------------------------------------
  observeEvent(input$nav_to, {
    updateTabsetPanel(session, "main_tabs", selected = input$nav_to)
  })

  ##----------------------------------------------------------------------
  ##  0. Put your data in a named list for easy reference
  ##----------------------------------------------------------------------
  dataset_list <- list(
    control  = control,
    lowinput = lowinput
  )

  ##----------------------------------------------------------------------
  ##  Helper function: Get dataset based on data type selection
  ##----------------------------------------------------------------------
  get_dataset <- function(dataset_name, data_type) {
    if (data_type %in% c("tic", "clr")) {
      dataset_list_tic[[dataset_name]]
    } else {
      dataset_list[[dataset_name]]
    }
  }

  data_type_suffix <- function(data_type) {
    switch(data_type, tic = "_TIC", clr = "_CLR", "_raw")
  }

  data_type_label <- function(data_type) {
    switch(data_type, tic = "%TIC", clr = "CLT/CLR", "Raw")
  }

  coerce_numeric_table <- function(tbl) {
    as.data.frame(lapply(tbl, function(col) {
      if (is.numeric(col)) {
        return(as.numeric(col))
      }
      col_chr <- gsub(",", "", trimws(as.character(col)))
      suppressWarnings(as.numeric(col_chr))
    }), stringsAsFactors = FALSE, check.names = FALSE)
  }

  compute_pca_bundle <- function(df, mode, data_type, class_filter = "All",
                                 center = TRUE, scale = TRUE, pseudocount = 1e-6) {
    if (is.null(df) || ncol(df) <= 1 || nrow(df) < 2) {
      return(NULL)
    }

    if (is.null(class_filter) || length(class_filter) == 0 || is.na(class_filter[[1]])) {
      class_filter <- "All"
    } else {
      class_filter <- class_filter[[1]]
    }

    numeric_df <- df[, -1, drop = FALSE]
    numeric_df <- coerce_numeric_table(numeric_df)

    if (mode == "Lipids" && !is.null(class_filter) && class_filter != "All") {
      lipid_classes <- infer_lipid_class(colnames(numeric_df))
      names(lipid_classes) <- colnames(numeric_df)
      selected_cols <- names(lipid_classes)[lipid_classes == class_filter]
      numeric_df <- numeric_df[, selected_cols, drop = FALSE]
    }

    if (ncol(numeric_df) < 2) {
      return(NULL)
    }

    X <- as.matrix(numeric_df)
    storage.mode(X) <- "numeric"
    X[!is.finite(X)] <- 0

    if (data_type == "clr") {
      # Compositional PCA: closure -> pseudocount -> CLR -> z-score -> PCA
      X[X < 0] <- 0

      row_totals <- rowSums(X, na.rm = TRUE)
      row_totals[row_totals <= 0 | !is.finite(row_totals)] <- NA_real_
      P <- X / row_totals
      P[!is.finite(P)] <- 0

      P <- P + pseudocount
      P <- P / rowSums(P, na.rm = TRUE)

      L <- log(P)
      X_clr <- L - rowMeans(L, na.rm = TRUE)

      if (mode == "Samples") {
        V <- scale(X_clr, center = TRUE, scale = TRUE)
      } else {
        V <- t(X_clr)
        V <- t(scale(t(V), center = TRUE, scale = TRUE))
      }

      V[!is.finite(V)] <- 0
      pca_res <- prcomp(V, center = FALSE, scale. = FALSE)
      data_label <- "CLT/CLR"
    } else {
      if (mode == "Samples") {
        pca_res <- prcomp(X, center = center, scale. = scale)
      } else {
        V <- t(X)
        colnames(V) <- df[[1]]
        keep_cols <- apply(V, 2, function(col) {
          finite_col <- col[is.finite(col)]
          length(unique(finite_col)) > 1
        })
        if (sum(keep_cols) < 2) {
          return(NULL)
        }
        V <- V[, keep_cols, drop = FALSE]
        pca_res <- prcomp(V, center = center, scale. = scale)
      }
      data_label <- data_type_label(data_type)
    }

    if (is.null(pca_res$x) || ncol(pca_res$x) < 2) {
      return(NULL)
    }

    pca_scores <- as.data.frame(pca_res$x[, 1:2, drop = FALSE])
    colnames(pca_scores) <- c("PC1", "PC2")

    if (mode == "Samples") {
      pca_scores$ID <- df[[1]]
    } else {
      pca_scores$ID <- rownames(pca_scores)
      pca_scores$Class <- infer_lipid_class(pca_scores$ID)
      pca_scores$Group <- class_to_group(pca_scores$Class)
    }

    pve <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
    list(
      scores = pca_scores,
      pve = pve,
      data_label = data_label
    )
  }

  build_class_palette <- function(classes) {
    base <- class_colors
    missing_classes <- setdiff(unique(classes), names(base))
    if (length(missing_classes) == 0) {
      return(base)
    }
    extra_cols <- grDevices::hcl.colors(length(missing_classes), palette = "Dark 3")
    c(base, setNames(extra_cols, missing_classes))
  }

  build_lipid_pca_plot <- function(pca_data, title_text, xlab, ylab) {
    class_palette <- build_class_palette(pca_data$Class)
    ellipse_df <- pca_data %>%
      filter(!is.na(Group)) %>%
      group_by(Group) %>%
      filter(n() >= 3) %>%
      ungroup()

    p <- ggplot(pca_data, aes(x = PC1, y = PC2))

    if (nrow(ellipse_df) > 0) {
      p <- p +
        stat_ellipse(
          data = ellipse_df,
          aes(color = Group, group = Group),
          level = 0.95,
          linewidth = 1.2,
          type = "norm",
          show.legend = FALSE
        ) +
        scale_color_manual(values = group_colors, na.value = "grey70") +
        ggnewscale::new_scale_color()
    }

    p +
      geom_point(aes(color = Class, text = ID), size = 3, alpha = 0.9) +
      scale_color_manual(values = class_palette, name = "Class") +
      labs(
        title = title_text,
        x = xlab,
        y = ylab
      ) +
      plot_theme +
      theme(legend.position = "right")
  }

  normalize_threshold_choice <- function(p_threshold) {
    if (is.null(p_threshold) || length(p_threshold) == 0 || is.na(p_threshold[[1]])) {
      p_chr <- "7"
    } else {
      p_chr <- as.character(p_threshold[[1]])
    }
    if (!p_chr %in% names(annotation_data_by_threshold)) {
      p_chr <- "7"
    }
    if (!p_chr %in% names(annotation_data_by_threshold)) {
      p_chr <- names(annotation_data_by_threshold)[[1]]
    }
    p_chr
  }

  get_annotation_data_list <- function(dataset_name, p_threshold = "7") {
    p_key <- normalize_threshold_choice(p_threshold)
    data_pool <- annotation_data_by_threshold[[p_key]]
    if (is.null(data_pool) || is.null(data_pool[[dataset_name]])) {
      return(list())
    }
    data_pool[[dataset_name]]
  }

  get_annotation_source_map <- function(dataset_name, p_threshold = "7") {
    p_key <- normalize_threshold_choice(p_threshold)
    source_pool <- annotation_data_source_by_threshold[[p_key]]
    if (is.null(source_pool) || is.null(source_pool[[dataset_name]])) {
      return(setNames(character(0), character(0)))
    }
    source_pool[[dataset_name]]
  }

  get_annotation_traits <- function(dataset_name, source_choice = "all", p_threshold = "7") {
    data_list <- get_annotation_data_list(dataset_name, p_threshold = p_threshold)
    if (is.null(data_list) || length(data_list) == 0) {
      return(character(0))
    }

    trait_names <- names(data_list)
    source_map <- get_annotation_source_map(dataset_name, p_threshold = p_threshold)

    if (!is.null(source_map) && source_choice != "all") {
      trait_names <- trait_names[source_map[trait_names] == source_choice]
    }

    trait_names
  }

  extract_trait_tokens <- function(trait_name) {
    trait_name <- as.character(trait_name)
    if (length(trait_name) == 0 || is.na(trait_name) || trait_name == "") {
      return(character(0))
    }

    tokens <- character(0)

    lead_token <- stringr::str_match(trait_name, "^([A-Za-z0-9]+)(?=\\()")[, 2]
    if (!is.na(lead_token) && lead_token != "") {
      tokens <- c(tokens, lead_token)
    }

    mod_token <- stringr::str_match(trait_name, "^([A-Za-z0-9]+)_mod_sub")[, 2]
    if (!is.na(mod_token) && mod_token != "") {
      tokens <- c(tokens, mod_token)
    }

    ratio_tokens <- regmatches(
      trait_name,
      gregexpr("(?<=^Sum_)[A-Za-z0-9]+|(?<=_over_)[A-Za-z0-9]+", trait_name, perl = TRUE)
    )[[1]]
    if (length(ratio_tokens) > 0 && !identical(ratio_tokens, character(0))) {
      tokens <- c(tokens, ratio_tokens)
    }

    tokens <- unique(tokens[!is.na(tokens) & tokens != ""])
    tokens
  }

  build_trait_token_lookup <- function(dataset_name) {
    meta <- gwas_lipid_class_map[[dataset_name]]
    if (is.null(meta) || nrow(meta) == 0) {
      return(list(classes = setNames(list(), character(0)), subclasses = setNames(list(), character(0))))
    }

    token_tbl <- meta %>%
      transmute(
        Token = stringr::str_extract(Original_Lipid, "^[A-Za-z0-9]+(?=\\()"),
        Class = ifelse(is.na(Class), "", trimws(Class)),
        Subclass = ifelse(is.na(Subclass), "", trimws(Subclass))
      ) %>%
      filter(!is.na(Token), Token != "")

    if (nrow(token_tbl) == 0) {
      return(list(classes = setNames(list(), character(0)), subclasses = setNames(list(), character(0))))
    }

    class_tbl <- token_tbl %>% filter(Class != "")
    subclass_tbl <- token_tbl %>% filter(Subclass != "")

    classes_lookup <- if (nrow(class_tbl) == 0) {
      setNames(list(), character(0))
    } else {
      split(class_tbl$Class, class_tbl$Token) %>%
        lapply(function(x) unique(as.character(x)))
    }

    subclasses_lookup <- if (nrow(subclass_tbl) == 0) {
      setNames(list(), character(0))
    } else {
      split(subclass_tbl$Subclass, subclass_tbl$Token) %>%
        lapply(function(x) unique(as.character(x)))
    }

    list(classes = classes_lookup, subclasses = subclasses_lookup)
  }

  get_gwas_trait_meta_table <- function(dataset_name, trait_names = NULL, source_choice = "all", p_threshold = "7") {
    data_list <- get_annotation_data_list(dataset_name, p_threshold = p_threshold)
    all_traits <- names(data_list)

    if (is.null(trait_names)) {
      trait_names <- all_traits
    }
    trait_names <- as.character(trait_names)
    if (length(trait_names) == 0) {
      return(tibble::tibble(
        Trait = character(),
        Source = character(),
        Class = character(),
        Subclass = character(),
        ClassMembers = list(),
        SubclassMembers = list()
      ))
    }

    source_map <- get_annotation_source_map(dataset_name, p_threshold = p_threshold)
    source_values <- if (is.null(source_map)) {
      setNames(rep("individual", length(trait_names)), trait_names)
    } else {
      source_map[trait_names]
    }
    source_values[is.na(source_values)] <- "unknown"

    trait_tbl <- tibble::tibble(
      Trait = trait_names,
      Source = as.character(source_values),
      Key = normalize_lipid_key(trait_names)
    )

    if (!is.null(source_choice) && source_choice != "all") {
      trait_tbl <- trait_tbl %>% filter(Source == source_choice)
      if (nrow(trait_tbl) == 0) {
        return(tibble::tibble(
          Trait = character(),
          Source = character(),
          Class = character(),
          Subclass = character(),
          ClassMembers = list(),
          SubclassMembers = list()
        ))
      }
    }

    meta <- gwas_lipid_class_map[[dataset_name]]
    direct_map <- if (is.null(meta) || nrow(meta) == 0) {
      tibble::tibble(Key = character(), DirectClass = list(), DirectSubclass = list())
    } else {
      meta %>%
        mutate(
          Key = normalize_lipid_key(Original_Lipid),
          Class = ifelse(is.na(Class), "", trimws(Class)),
          Subclass = ifelse(is.na(Subclass), "", trimws(Subclass))
        ) %>%
        filter(Key != "") %>%
        group_by(Key) %>%
        summarise(
          DirectClass = list(unique(Class[Class != ""])),
          DirectSubclass = list(unique(Subclass[Subclass != ""])),
          .groups = "drop"
        )
    }

    trait_tbl <- trait_tbl %>%
      left_join(direct_map, by = "Key")

    trait_tbl$DirectClass <- lapply(trait_tbl$DirectClass, function(x) {
      if (is.null(x) || (length(x) == 1 && is.na(x))) character(0) else as.character(x)
    })
    trait_tbl$DirectSubclass <- lapply(trait_tbl$DirectSubclass, function(x) {
      if (is.null(x) || (length(x) == 1 && is.na(x))) character(0) else as.character(x)
    })

    token_lookup <- build_trait_token_lookup(dataset_name)

    class_members <- mapply(
      function(trait, direct_class) {
        tokens <- extract_trait_tokens(trait)
        token_classes <- unique(unlist(token_lookup$classes[tokens], use.names = FALSE))
        out <- unique(c(direct_class, token_classes))
        out[!is.na(out) & out != ""]
      },
      trait_tbl$Trait,
      trait_tbl$DirectClass,
      SIMPLIFY = FALSE
    )

    subclass_members <- mapply(
      function(trait, direct_subclass) {
        tokens <- extract_trait_tokens(trait)
        token_subclasses <- unique(unlist(token_lookup$subclasses[tokens], use.names = FALSE))
        out <- unique(c(direct_subclass, token_subclasses))
        out[!is.na(out) & out != ""]
      },
      trait_tbl$Trait,
      trait_tbl$DirectSubclass,
      SIMPLIFY = FALSE
    )

    class_display <- vapply(class_members, function(x) {
      if (length(x) == 0) "" else paste(x, collapse = " | ")
    }, character(1))

    subclass_display <- vapply(subclass_members, function(x) {
      if (length(x) == 0) "" else if (length(x) == 1) x[[1]] else "Multiple"
    }, character(1))

    tibble::tibble(
      Trait = trait_tbl$Trait,
      Source = trait_tbl$Source,
      Class = class_display,
      Subclass = subclass_display,
      ClassMembers = class_members,
      SubclassMembers = subclass_members
    )
  }

  get_filtered_gwas_traits <- function(dataset_name, source_choice = "all", class_choice = "All", subclass_choice = "All", lipid_choice = "All", p_threshold = "7") {
    meta_tbl <- get_gwas_trait_meta_table(dataset_name, source_choice = source_choice, p_threshold = p_threshold)
    if (nrow(meta_tbl) == 0) {
      return(character(0))
    }

    if (!is.null(class_choice) && class_choice != "All") {
      meta_tbl <- meta_tbl[vapply(meta_tbl$ClassMembers, function(x) class_choice %in% x, logical(1)), , drop = FALSE]
    }

    if (!is.null(subclass_choice) && subclass_choice != "All") {
      meta_tbl <- meta_tbl[vapply(meta_tbl$SubclassMembers, function(x) subclass_choice %in% x, logical(1)), , drop = FALSE]
    }

    if (!is.null(lipid_choice) && lipid_choice != "All") {
      meta_tbl <- meta_tbl %>% filter(Trait == lipid_choice)
    }

    unique(meta_tbl$Trait)
  }

  ##======================================================================
  ##  DOWNLOAD HANDLERS
  ##======================================================================

  # Reactive values to store plots for downloading
  plot_store <- reactiveValues(
    boxplot_single = NULL,
    boxplot_both = NULL,
    boxplot_ratio = NULL,
    boxplot_ratio_combined = NULL,
    correlation_plot = NULL,
    geno_view_plot = NULL,
    venn_plot = NULL
  )

  # Data Preview - Download Table

  output$download_data_table <- downloadHandler(
    filename = function() {
      class_suffix <- if (!is.null(input$class_choice) && input$class_choice != "All") {
        paste0("_", input$class_choice)
      } else { "" }
      subclass_suffix <- if (!is.null(input$subclass_choice_data) && input$subclass_choice_data != "All") {
        paste0("_", gsub("[^A-Za-z0-9]+", "_", input$subclass_choice_data))
      } else { "" }
      lipid_suffix <- if (!is.null(input$lipid_choice_data) && input$lipid_choice_data != "All") {
        paste0("_", gsub("[^A-Za-z0-9]+", "_", input$lipid_choice_data))
      } else { "" }
      type_suffix <- data_type_suffix(input$data_type_data)
      paste0(input$data_choice, "_lipids", class_suffix, subclass_suffix, lipid_suffix, type_suffix, ".csv")
    },
    content = function(file) {
      df <- get_dataset(input$data_choice, input$data_type_data)
      df_t <- as.data.frame(t(df[, -1]))
      colnames(df_t) <- df[[1]]
      df_t$LipidName <- rownames(df_t)
      rownames(df_t) <- NULL
      df_t <- df_t[, c(ncol(df_t), 1:(ncol(df_t)-1))]
      allowed_lipids <- get_data_preview_allowed_lipids()
      if (length(allowed_lipids) > 0) {
        df_t <- df_t %>% filter(LipidName %in% allowed_lipids)
      }
      write.csv(df_t, file, row.names = FALSE)
    }
  )

  # GWAS - Download Table
  output$download_gwas_table <- downloadHandler(
    filename = function() {
      paste0("GWAS_", input$gwas_dataset, "_plog10", input$gwas_threshold, "_", input$selected_lipid, ".csv")
    },
    content = function(file) {
      data_list <- get_annotation_data_list(input$gwas_dataset, p_threshold = input$gwas_threshold)
      if (input$selected_lipid == "__ALL__") {
        traits <- get_filtered_gwas_traits(
          input$gwas_dataset,
          source_choice = input$gwas_trait_source,
          class_choice = input$gwas_class_choice,
          subclass_choice = input$gwas_subclass_choice,
          p_threshold = input$gwas_threshold
        )
        rows <- lapply(traits, function(trait) {
          df <- data_list[[trait]]
          if (is.null(df) || nrow(df) == 0) {
            return(NULL)
          }
          meta <- get_gwas_trait_meta_table(
            input$gwas_dataset,
            trait_names = trait,
            source_choice = input$gwas_trait_source,
            p_threshold = input$gwas_threshold
          )
          df %>% mutate(
            Trait = trait,
            Class = meta$Class[[1]],
            Subclass = meta$Subclass[[1]],
            .before = GeneID
          )
        })
        rows <- rows[!vapply(rows, is.null, logical(1))]
        if (length(rows) == 0) {
          write.csv(data.frame(), file, row.names = FALSE)
          return()
        }
        df_annotated <- bind_rows(rows) %>% left_join(gene_annotation, by = "GeneID")
        write.csv(df_annotated, file, row.names = FALSE)
        return()
      }

      df <- data_list[[input$selected_lipid]]
      if (is.null(df)) {
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }

      df_annotated <- df %>% left_join(gene_annotation, by = "GeneID")
      meta <- get_gwas_trait_meta_table(
        input$gwas_dataset,
        trait_names = input$selected_lipid,
        source_choice = input$gwas_trait_source,
        p_threshold = input$gwas_threshold
      )
      if (!is.null(meta) && nrow(meta) > 0) {
        df_annotated <- df_annotated %>%
          mutate(
            Class = meta$Class[[1]],
            Subclass = meta$Subclass[[1]],
            .before = GeneID
          )
      }
      write.csv(df_annotated, file, row.names = FALSE)
    }
  )

  # Gene Hits - Download Table
  output$download_gene_hits <- downloadHandler(
    filename = function() {
      class_suffix <- if (!is.null(input$hit_class_choice) && input$hit_class_choice != "All") {
        paste0("_", input$hit_class_choice)
      } else { "" }
      paste0("gene_hits_", input$hit_dataset, "_plog10", input$hit_threshold, class_suffix, ".csv")
    },
    content = function(file) {
      data_list <- get_annotation_data_list(input$hit_dataset, p_threshold = input$hit_threshold)
      if (length(data_list) == 0) {
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }
      threshold <- as.numeric(input$hit_threshold)
      trait_subset <- get_filtered_gwas_traits(
        input$hit_dataset,
        source_choice = input$hit_trait_source,
        class_choice = input$hit_class_choice,
        subclass_choice = input$hit_subclass_choice,
        lipid_choice = input$hit_lipid_choice,
        p_threshold = input$hit_threshold
      )
      hit_rows <- lapply(trait_subset, function(trait) {
        df <- data_list[[trait]]
        if (!is.null(df) && all(c("GeneID", "log(p)") %in% colnames(df))) {
          out <- df[df$`log(p)` >= threshold, c("GeneID", "log(p)"), drop = FALSE]
          out$Trait <- trait
          out
        }
      })
      hit_rows <- hit_rows[!vapply(hit_rows, is.null, logical(1))]
      if (length(hit_rows) == 0) {
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }
      combined_df <- dplyr::bind_rows(hit_rows)
      gene_summary <- combined_df %>%
        group_by(GeneID) %>%
        summarise(
          Occurrences = n(),
          Highest_logp = round(max(`log(p)`, na.rm = TRUE), 2),
          .groups = "drop"
        ) %>%
        arrange(desc(Occurrences), desc(Highest_logp)) %>%
        left_join(gene_annotation, by = "GeneID")
      write.csv(gene_summary, file, row.names = FALSE)
    }
  )

  # Phenotype Hits - Download Table
  output$download_phenotype_hits <- downloadHandler(
    filename = function() {
      class_suffix <- if (!is.null(input$hit_class_choice) && input$hit_class_choice != "All") {
        paste0("_", input$hit_class_choice)
      } else { "" }
      paste0("phenotype_hits_", input$hit_dataset, "_plog10", input$hit_threshold, class_suffix, ".csv")
    },
    content = function(file) {
      data_list <- get_annotation_data_list(input$hit_dataset, p_threshold = input$hit_threshold)
      if (length(data_list) == 0) {
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }

      threshold <- as.numeric(input$hit_threshold)
      trait_subset <- get_filtered_gwas_traits(
        input$hit_dataset,
        source_choice = input$hit_trait_source,
        class_choice = input$hit_class_choice,
        subclass_choice = input$hit_subclass_choice,
        lipid_choice = input$hit_lipid_choice,
        p_threshold = input$hit_threshold
      )
      trait_meta <- get_gwas_trait_meta_table(
        input$hit_dataset,
        trait_names = trait_subset,
        source_choice = input$hit_trait_source,
        p_threshold = input$hit_threshold
      ) %>%
        transmute(Trait = Trait, Class = Class, Subclass = Subclass)
      source_map <- get_annotation_source_map(input$hit_dataset, p_threshold = input$hit_threshold)

      summary_rows <- lapply(trait_subset, function(trait) {
        df <- data_list[[trait]]
        if (is.null(df) || !all(c("GeneID", "log(p)") %in% colnames(df))) {
          return(NULL)
        }
        df_sig <- df[df$`log(p)` >= threshold, c("GeneID", "log(p)"), drop = FALSE]
        if (nrow(df_sig) == 0) {
          return(NULL)
        }
        idx <- which.max(df_sig$`log(p)`)
        trait_info <- trait_meta %>% filter(Trait == trait)
        data.frame(
          Phenotype = trait,
          TraitSource = if (!is.null(source_map)) source_map[[trait]] else "unknown",
          Class = if (nrow(trait_info) > 0) trait_info$Class[[1]] else "",
          Subclass = if (nrow(trait_info) > 0) trait_info$Subclass[[1]] else "",
          SignificantGenes = nrow(df_sig),
          Highest_logp = round(max(df_sig$`log(p)`, na.rm = TRUE), 2),
          TopGeneID = df_sig$GeneID[[idx]],
          stringsAsFactors = FALSE
        )
      })
      summary_rows <- summary_rows[!vapply(summary_rows, is.null, logical(1))]

      if (length(summary_rows) == 0) {
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }

      out <- bind_rows(summary_rows) %>%
        arrange(desc(Highest_logp), desc(SignificantGenes), Phenotype)
      write.csv(out, file, row.names = FALSE)
    }
  )

  # Lipid Lookup Table - Download
  output$download_lipid_lookup <- downloadHandler(
    filename = function() {
      paste0("lipid_key_", input$heatmap_dataset, ".csv")
    },
    content = function(file) {
      df <- dataset_list[[input$heatmap_dataset]]
      mat <- df[, -1]
      original_lipid_names <- colnames(mat)
      short_lipid_ids <- paste0("L", seq_along(original_lipid_names))
      lipid_name_map <- data.frame(ID = short_lipid_ids, Name = original_lipid_names)
      write.csv(lipid_name_map, file, row.names = FALSE)
    }
  )

  # Lipid Lookup Table for Condition Comparison - Download
  output$download_lipid_lookup_diff <- downloadHandler(
    filename = function() {
      paste0("lipid_key_condition_comparison.csv")
    },
    content = function(file) {
      df_control <- get_dataset("control", input$data_type_heatmap_diff)
      df_lowinput <- get_dataset("lowinput", input$data_type_heatmap_diff)
      all_lipids <- sort(union(colnames(df_control)[-1], colnames(df_lowinput)[-1]))
      short_lipid_ids <- paste0("L", seq_along(all_lipids))
      lipid_name_map <- data.frame(ID = short_lipid_ids, Name = all_lipids)
      write.csv(lipid_name_map, file, row.names = FALSE)
    }
  )

  # Venn Table - Download
  output$download_venn_table <- downloadHandler(
    filename = function() {
      type_suffix <- if (input$venn_type == "Class" && !is.null(input$venn_class_choice)) {
        paste0("_", input$venn_class_choice)
      } else { "_all" }
      paste0("venn_sets", type_suffix, ".csv")
    },
    content = function(file) {
      df_control <- dataset_list$control
      df_lowinput <- dataset_list$lowinput
      if (input$venn_type == "Individual") {
        set_control <- colnames(df_control)[-1]
        set_lowinput <- colnames(df_lowinput)[-1]
      } else {
        class_lipids <- lipid_class_info$Lipids[lipid_class_info$Class == input$venn_class_choice]
        set_control <- intersect(class_lipids, colnames(df_control))
        set_lowinput <- intersect(class_lipids, colnames(df_lowinput))
      }
      only_control <- setdiff(set_control, set_lowinput)
      only_lowinput <- setdiff(set_lowinput, set_control)
      shared <- intersect(set_control, set_lowinput)
      max_len <- max(length(only_control), length(shared), length(only_lowinput))
      venn_df <- data.frame(
        Only_Control = c(only_control, rep(NA, max_len - length(only_control))),
        Shared = c(shared, rep(NA, max_len - length(shared))),
        Only_Low_Input = c(only_lowinput, rep(NA, max_len - length(only_lowinput)))
      )
      write.csv(venn_df, file, row.names = FALSE)
    }
  )

  # ---- PLOT DOWNLOAD HANDLERS ----

  # Boxplot Single (Density Plot)
  output$download_boxplot_single <- downloadHandler(
    filename = function() {
      compound <- if (input$compound_choice == "ALL_LIPIDS_SUM") input$class_filter_plot else input$compound_choice
      paste0("density_", input$plot_dataset, "_", compound, ".png")
    },
    content = function(file) {
      p <- plot_store$boxplot_single
      if (!is.null(p)) ggsave(file, p, width = 10, height = 6, dpi = 300)
    }
  )

  # Boxplot Both (Comparison Plot)
  output$download_boxplot_both <- downloadHandler(
    filename = function() {
      compound <- if (input$compound_choice == "ALL_LIPIDS_SUM") input$class_filter_plot else input$compound_choice
      paste0("comparison_", compound, ".png")
    },
    content = function(file) {
      p <- plot_store$boxplot_both
      if (!is.null(p)) ggsave(file, p, width = 10, height = 6, dpi = 300)
    }
  )

  # Ratio Single
  output$download_ratio_single <- downloadHandler(
    filename = function() {
      paste0("ratio_", input$ratio_lipid_1, "_", input$ratio_lipid_2, "_", input$plot_dataset_ratio, ".png")
    },
    content = function(file) {
      p <- plot_store$boxplot_ratio
      if (!is.null(p)) ggsave(file, p, width = 10, height = 6, dpi = 300)
    }
  )

  # Ratio Combined
  output$download_ratio_combined <- downloadHandler(
    filename = function() {
      paste0("ratio_combined_", input$ratio_lipid_1, "_", input$ratio_lipid_2, ".png")
    },
    content = function(file) {
      p <- plot_store$boxplot_ratio_combined
      if (!is.null(p)) ggsave(file, p, width = 10, height = 6, dpi = 300)
    }
  )

  # PCA Plot (using plotly export)
  output$download_pca_plot <- downloadHandler(
    filename = function() {
      type_suffix <- data_type_suffix(input$data_type_pca)
      paste0("PCA_", input$pca_dataset, "_", input$pca_mode, type_suffix, ".png")
    },
    content = function(file) {
      pca_data_mode <- if (input$data_type_pca == "clr") "tic" else input$data_type_pca
      df <- get_dataset(input$pca_dataset, pca_data_mode)
      pca_bundle <- compute_pca_bundle(
        df = df,
        mode = input$pca_mode,
        data_type = input$data_type_pca,
        class_filter = if (input$pca_mode == "Lipids") input$pca_class_filter else "All",
        center = input$pca_center,
        scale = input$pca_scale
      )

      if (is.null(pca_bundle)) {
        p <- ggplot() +
          theme_void() +
          annotate("text", x = 0.5, y = 0.5, label = "PCA not available for this selection.")
      } else {
        pca_data <- pca_bundle$scores
        xlab <- sprintf("PC1 (%.1f%%)", 100 * pca_bundle$pve[1])
        ylab <- sprintf("PC2 (%.1f%%)", 100 * pca_bundle$pve[2])

        if (input$pca_mode == "Samples") {
          p <- ggplot(pca_data, aes(x = PC1, y = PC2, label = ID)) +
            geom_point(size = 3, alpha = 0.7, color = "#667eea") +
            labs(
              title = paste("PCA - Samples (", input$pca_dataset, ",", pca_bundle$data_label, ")"),
              x = xlab,
              y = ylab
            ) +
            plot_theme
        } else {
          p <- build_lipid_pca_plot(
            pca_data = pca_data,
            title_text = paste("PCA - Lipids (", input$pca_dataset, ",", pca_bundle$data_label, ")"),
            xlab = xlab,
            ylab = ylab
          )
        }
      }
      ggsave(file, p, width = 12, height = 8, dpi = 300)
    }
  )

  # t-SNE/UMAP Plot
  output$download_dr_plot <- downloadHandler(
    filename = function() {
      type_suffix <- data_type_suffix(input$data_type_dr)
      paste0(input$dr_method, "_", input$dr_dataset, "_", input$dr_mode, type_suffix, ".png")
    },
    content = function(file) {
      df <- get_dataset(input$dr_dataset, input$data_type_dr)
      method <- input$dr_method
      data_label <- data_type_label(input$data_type_dr)

      if (input$dr_mode == "Samples") {
        numeric_df <- df[, -1]
        numeric_df <- as.data.frame(sapply(numeric_df, as.numeric))
        row_labels <- df[[1]]

        if (method == "t-SNE") {
          dr_out <- Rtsne(numeric_df, perplexity = input$perplexity, check_duplicates = FALSE)
          dim_df <- data.frame(Dim1 = dr_out$Y[,1], Dim2 = dr_out$Y[,2], Label = row_labels)
        } else {
          dr_out <- umap(numeric_df)
          dim_df <- data.frame(Dim1 = dr_out$layout[,1], Dim2 = dr_out$layout[,2], Label = row_labels)
        }

        p <- ggplot(dim_df, aes(x = Dim1, y = Dim2)) +
          geom_point(size = 3, alpha = 0.7, color = "#667eea") +
          labs(title = paste(method, "- Samples (", input$dr_dataset, ",", data_label, ")"), x = "Dim 1", y = "Dim 2") +
          plot_theme

      } else {
        # Lipids mode with ellipses
        numeric_df <- df[, -1]
        numeric_df <- t(numeric_df)
        colnames(numeric_df) <- df[[1]]
        numeric_df <- as.data.frame(numeric_df)
        row_labels <- rownames(numeric_df)

        # Extract class from lipid name using get_lipid_class (like paper script)
        class_col <- get_lipid_class(row_labels)

        if (method == "t-SNE") {
          dr_out <- Rtsne(numeric_df, perplexity = input$perplexity, check_duplicates = FALSE)
          dim_df <- data.frame(Dim1 = dr_out$Y[,1], Dim2 = dr_out$Y[,2], Label = row_labels, Class = class_col)
        } else {
          dr_out <- umap(numeric_df)
          dim_df <- data.frame(Dim1 = dr_out$layout[,1], Dim2 = dr_out$layout[,2], Label = row_labels, Class = class_col)
        }

        dim_df$Group <- class_to_group(dim_df$Class)

        # Filter out "Other" class lipids
        dim_df <- dim_df[dim_df$Class != "Other", ]
        dim_df <- dim_df[!is.na(dim_df$Group), ]

        # ggplot with ellipses (exactly like paper script)
        p <- ggplot(dim_df, aes(x = Dim1, y = Dim2)) +
          stat_ellipse(aes(color = Group), level = 0.95, linewidth = 1.5, type = "norm", show.legend = FALSE) +
          scale_color_manual(values = group_colors) +
          ggnewscale::new_scale_color() +
          geom_point(aes(color = Class), size = 3, alpha = 0.9) +
          scale_color_manual(values = class_colors, name = "Class") +
          labs(title = paste(method, "- Lipids (", input$dr_dataset, ",", data_label, ")"), x = "Dim 1", y = "Dim 2") +
          plot_theme +
          theme(legend.position = "right")
      }
      ggsave(file, p, width = 12, height = 8, dpi = 300)
    }
  )

  # Correlation Plot
  output$download_correlation_plot <- downloadHandler(
    filename = function() {
      paste0("correlation_", input$cor_dataset, "_", input$lipid1, "_vs_", input$lipid2, ".png")
    },
    content = function(file) {
      p <- plot_store$correlation_plot
      if (!is.null(p)) ggsave(file, p, width = 10, height = 8, dpi = 300)
    }
  )

  # Genotype View Plot
  output$download_geno_plot <- downloadHandler(
    filename = function() {
      paste0("genotype_", input$geno_view_id, "_", input$geno_view_lipid, ".png")
    },
    content = function(file) {
      p <- plot_store$geno_view_plot
      if (!is.null(p)) ggsave(file, p, width = 8, height = 6, dpi = 300)
    }
  )

  # Venn Plot
  output$download_venn_plot <- downloadHandler(
    filename = function() {
      type_suffix <- if (input$venn_type == "Class" && !is.null(input$venn_class_choice)) {
        paste0("_", input$venn_class_choice)
      } else { "_all" }
      paste0("venn_diagram", type_suffix, ".png")
    },
    content = function(file) {
      df_control <- dataset_list$control
      df_lowinput <- dataset_list$lowinput
      if (input$venn_type == "Individual") {
        set_control <- colnames(df_control)[-1]
        set_lowinput <- colnames(df_lowinput)[-1]
      } else {
        class_lipids <- lipid_class_info$Lipids[lipid_class_info$Class == input$venn_class_choice]
        set_control <- intersect(class_lipids, colnames(df_control))
        set_lowinput <- intersect(class_lipids, colnames(df_lowinput))
      }
      png(file, width = 800, height = 800, res = 150)
      venn.plot <- draw.pairwise.venn(
        area1 = length(set_control),
        area2 = length(set_lowinput),
        cross.area = length(intersect(set_control, set_lowinput)),
        category = c("Control", "Low Input"),
        fill = c("skyblue", "orange"),
        alpha = c(0.5, 0.5),
        lty = "blank",
        cex = 2,
        cat.cex = 1.5,
        cat.pos = c(-20, 20),
        cat.dist = 0.05,
        scaled = TRUE
      )
      grid.draw(venn.plot)
      dev.off()
    }
  )

  ##======================================================================
  ##  TAB 1: Data Preview
  ##======================================================================

  data_subclass_col <- dplyr::case_when(
    "Subclass" %in% names(lipid_class_info) ~ "Subclass",
    "SubClass" %in% names(lipid_class_info) ~ "SubClass",
    "Secondary_Class" %in% names(lipid_class_info) ~ "Secondary_Class",
    TRUE ~ NA_character_
  )

  get_data_preview_meta <- function(dataset_name) {
    df <- dataset_list[[dataset_name]]
    if (is.null(df) || ncol(df) <= 1) {
      return(tibble::tibble(LipidName = character(), Class = character(), Subclass = character()))
    }

    lipids_in_data <- colnames(df)[-1]
    base_tbl <- tibble::tibble(
      LipidName = lipids_in_data,
      Key = normalize_lipid_key(lipids_in_data)
    )

    if (is.null(lipid_class_info) || nrow(lipid_class_info) == 0 || !"Lipids" %in% names(lipid_class_info)) {
      return(base_tbl %>% mutate(Class = "Unclassified", Subclass = "Unclassified") %>% select(LipidName, Class, Subclass))
    }

    map_tbl <- lipid_class_info %>%
      transmute(
        LipidName = as.character(Lipids),
        Class = as.character(Class),
        Subclass = if (!is.na(data_subclass_col)) as.character(.data[[data_subclass_col]]) else ""
      ) %>%
      mutate(
        Key = normalize_lipid_key(LipidName),
        Class = ifelse(is.na(Class) | trimws(Class) == "", "Unclassified", trimws(Class)),
        Subclass = ifelse(is.na(Subclass) | trimws(Subclass) == "", "Unclassified", trimws(Subclass))
      ) %>%
      filter(Key != "") %>%
      distinct(Key, .keep_all = TRUE)

    base_tbl %>%
      left_join(map_tbl %>% select(Key, Class, Subclass), by = "Key") %>%
      mutate(
        Class = ifelse(is.na(Class) | Class == "", "Unclassified", Class),
        Subclass = ifelse(is.na(Subclass) | Subclass == "", "Unclassified", Subclass)
      ) %>%
      select(LipidName, Class, Subclass)
  }

  get_data_preview_allowed_lipids <- reactive({
    req(input$data_choice)
    meta <- get_data_preview_meta(input$data_choice)

    if (!is.null(input$class_choice) && input$class_choice != "All") {
      meta <- meta %>% filter(Class == input$class_choice)
    }
    if (!is.null(input$subclass_choice_data) && input$subclass_choice_data != "All") {
      meta <- meta %>% filter(Subclass == input$subclass_choice_data)
    }
    if (!is.null(input$lipid_choice_data) && input$lipid_choice_data != "All") {
      meta <- meta %>% filter(LipidName == input$lipid_choice_data)
    }

    unique(meta$LipidName)
  })
  
  output$class_filter_ui <- renderUI({
    req(input$data_choice)
    meta <- get_data_preview_meta(input$data_choice)
    classes <- sort(unique(meta$Class))
    selectInput(
      inputId = "class_choice",
      label = "Filter by Lipid Class:",
      choices = c("All", classes),
      selected = "All"
    )
  })

  output$subclass_filter_data_ui <- renderUI({
    req(input$data_choice, input$class_choice)
    meta <- get_data_preview_meta(input$data_choice)
    if (input$class_choice != "All") {
      meta <- meta %>% filter(Class == input$class_choice)
    }
    subclasses <- sort(unique(meta$Subclass))
    selectInput(
      inputId = "subclass_choice_data",
      label = "Filter by Subclass:",
      choices = c("All", subclasses),
      selected = "All"
    )
  })

  output$lipid_filter_data_ui <- renderUI({
    req(input$data_choice, input$class_choice, input$subclass_choice_data)
    meta <- get_data_preview_meta(input$data_choice)
    if (input$class_choice != "All") {
      meta <- meta %>% filter(Class == input$class_choice)
    }
    if (input$subclass_choice_data != "All") {
      meta <- meta %>% filter(Subclass == input$subclass_choice_data)
    }
    lipids <- sort(unique(meta$LipidName))
    selectInput(
      inputId = "lipid_choice_data",
      label = "Filter by Lipid:",
      choices = c("All", lipids),
      selected = "All"
    )
  })
  
  output$data_table <- renderDT({
    req(input$data_choice, input$data_type_data)
    tryCatch({
      df <- get_dataset(input$data_choice, input$data_type_data)
      validate(need(!is.null(df) && ncol(df) > 1, "Dataset is loaded but has no plottable lipid columns."))

      # Transpose for display
      df_t <- as.data.frame(t(df[, -1]))
      colnames(df_t) <- df[[1]]
      df_t$LipidName <- rownames(df_t)
      rownames(df_t) <- NULL
      df_t <- df_t[, c(ncol(df_t), 1:(ncol(df_t)-1))]

      allowed_lipids <- get_data_preview_allowed_lipids()
      if (length(allowed_lipids) > 0) {
        df_t <- df_t %>% filter(LipidName %in% allowed_lipids)
      }

      validate(need(nrow(df_t) > 0, "No lipids found for the selected class filter."))

      dt <- datatable(
        df_t,
        options = list(scrollX = TRUE, pageLength = 10, lengthMenu = c(5, 10, 25, 50)),
        rownames = FALSE
      )

      if (ncol(df_t) >= 2) {
        dt <- dt %>% formatRound(columns = 2:ncol(df_t), digits = 2)
      }

      dt
    }, error = function(e) {
      datatable(
        data.frame(Message = paste("Data table failed to render:", conditionMessage(e))),
        options = list(dom = "t"),
        rownames = FALSE
      )
    })
  })

  data_viz_selected_datasets <- reactive({
    req(input$data_viz_mode)
    if (input$data_viz_mode == "compare") {
      c("control", "lowinput")
    } else {
      input$data_viz_dataset
    }
  })

  data_viz_meta <- reactive({
    req(input$data_viz_domain, input$data_viz_mode)
    datasets <- data_viz_selected_datasets()

    if (input$data_viz_domain == "gwas") {
      req(input$data_viz_trait_source, input$data_viz_gwas_threshold)
      rows <- lapply(datasets, function(ds) {
        tbl <- get_gwas_trait_meta_table(
          ds,
          source_choice = input$data_viz_trait_source,
          p_threshold = input$data_viz_gwas_threshold
        )
        if (nrow(tbl) == 0) {
          return(NULL)
        }
        tbl %>%
          transmute(
            Dataset = ds,
            Trait = Trait,
            Class = ifelse(Class == "", "Unclassified", Class),
            Subclass = ifelse(Subclass == "", "Unclassified", Subclass),
            ClassMembers = lapply(ClassMembers, function(x) if (length(x) == 0) "Unclassified" else x),
            SubclassMembers = lapply(SubclassMembers, function(x) if (length(x) == 0) "Unclassified" else x)
          )
      })
      rows <- rows[!vapply(rows, is.null, logical(1))]
      if (length(rows) == 0) {
        return(tibble::tibble(
          Dataset = character(),
          Trait = character(),
          Class = character(),
          Subclass = character(),
          ClassMembers = list(),
          SubclassMembers = list()
        ))
      }
      return(bind_rows(rows))
    }

    rows <- lapply(datasets, function(ds) {
      m <- get_data_preview_meta(ds)
      if (nrow(m) == 0) {
        return(NULL)
      }
      m %>%
        transmute(
          Dataset = ds,
          Trait = LipidName,
          Class = ifelse(Class == "", "Unclassified", Class),
          Subclass = ifelse(Subclass == "", "Unclassified", Subclass),
          ClassMembers = lapply(Class, function(x) x),
          SubclassMembers = lapply(Subclass, function(x) x)
        )
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) {
      return(tibble::tibble(
        Dataset = character(),
        Trait = character(),
        Class = character(),
        Subclass = character(),
        ClassMembers = list(),
        SubclassMembers = list()
      ))
    }
    bind_rows(rows)
  })

  output$data_viz_class_filter_ui <- renderUI({
    meta <- data_viz_meta()
    classes <- sort(unique(unlist(meta$ClassMembers, use.names = FALSE)))
    classes <- classes[!is.na(classes) & classes != ""]
    selectInput("data_viz_class_choice", "Filter by Class:", choices = c("All", classes), selected = "All")
  })

  output$data_viz_subclass_filter_ui <- renderUI({
    req(input$data_viz_class_choice)
    meta <- data_viz_meta()
    if (input$data_viz_class_choice != "All") {
      keep <- vapply(meta$ClassMembers, function(x) input$data_viz_class_choice %in% x, logical(1))
      meta <- meta[keep, , drop = FALSE]
    }
    subclasses <- sort(unique(unlist(meta$SubclassMembers, use.names = FALSE)))
    subclasses <- subclasses[!is.na(subclasses) & subclasses != ""]
    selectInput("data_viz_subclass_choice", "Filter by Subclass:", choices = c("All", subclasses), selected = "All")
  })

  output$data_viz_lipid_filter_ui <- renderUI({
    req(input$data_viz_class_choice, input$data_viz_subclass_choice)
    meta <- data_viz_meta()
    if (input$data_viz_class_choice != "All") {
      keep <- vapply(meta$ClassMembers, function(x) input$data_viz_class_choice %in% x, logical(1))
      meta <- meta[keep, , drop = FALSE]
    }
    if (input$data_viz_subclass_choice != "All") {
      keep <- vapply(meta$SubclassMembers, function(x) input$data_viz_subclass_choice %in% x, logical(1))
      meta <- meta[keep, , drop = FALSE]
    }
    lipids <- sort(unique(meta$Trait))
    selectInput("data_viz_lipid_choice", "Filter by Lipid:", choices = c("All", lipids), selected = "All")
  })

  data_viz_filtered_meta <- reactive({
    req(input$data_viz_class_choice, input$data_viz_subclass_choice, input$data_viz_lipid_choice)
    meta <- data_viz_meta()
    if (nrow(meta) == 0) {
      return(meta)
    }

    if (input$data_viz_class_choice != "All") {
      keep <- vapply(meta$ClassMembers, function(x) input$data_viz_class_choice %in% x, logical(1))
      meta <- meta[keep, , drop = FALSE]
    }
    if (input$data_viz_subclass_choice != "All") {
      keep <- vapply(meta$SubclassMembers, function(x) input$data_viz_subclass_choice %in% x, logical(1))
      meta <- meta[keep, , drop = FALSE]
    }
    if (input$data_viz_lipid_choice != "All") {
      meta <- meta %>% filter(Trait == input$data_viz_lipid_choice)
    }
    meta
  })

  effective_data_viz_level <- reactive({
    req(input$data_viz_level, input$data_viz_class_choice, input$data_viz_subclass_choice, input$data_viz_lipid_choice)
    level <- input$data_viz_level

    # If class is fixed to one value and grouping is by class, drill down to subclasses.
    if (
      level == "Class" &&
      !is.null(input$data_viz_class_choice) && input$data_viz_class_choice != "All" &&
      !is.null(input$data_viz_subclass_choice) && input$data_viz_subclass_choice == "All" &&
      !is.null(input$data_viz_lipid_choice) && input$data_viz_lipid_choice == "All"
    ) {
      return("Subclass")
    }

    level
  })

  output$data_visualize_plot <- renderPlot({
    req(input$data_viz_domain, input$data_viz_level, input$data_viz_top_n)
    meta <- data_viz_filtered_meta()
    validate(need(nrow(meta) > 0, "No entries found for the selected filters."))

    datasets <- data_viz_selected_datasets()
    level <- effective_data_viz_level()

    summary_df <- if (input$data_viz_domain == "abundance") {
      req(input$data_viz_data_type)
      rows <- lapply(datasets, function(ds) {
        df <- get_dataset(ds, input$data_viz_data_type)
        if (is.null(df) || ncol(df) <= 1) {
          return(NULL)
        }

        num_df <- coerce_numeric_table(df[, -1, drop = FALSE])
        lipid_means <- colMeans(as.matrix(num_df), na.rm = TRUE)
        value_tbl <- tibble::tibble(Dataset = ds, Trait = names(lipid_means), Value = as.numeric(lipid_means))
        keep_tbl <- meta %>% filter(Dataset == ds) %>% select(Dataset, Trait, Class, Subclass)
        value_tbl %>% inner_join(keep_tbl, by = c("Dataset", "Trait"))
      })
      rows <- rows[!vapply(rows, is.null, logical(1))]
      validate(need(length(rows) > 0, "No numeric lipid matrix available for visualization."))
      joined <- bind_rows(rows)
      validate(need(nrow(joined) > 0, "No abundance values left after filtering."))

      if (level == "Class") {
        joined %>% mutate(Group = Class)
      } else if (level == "Subclass") {
        joined %>% mutate(Group = Subclass)
      } else {
        joined %>% mutate(Group = Trait)
      }
    } else {
      rows <- lapply(seq_len(nrow(meta)), function(i) {
        ds <- meta$Dataset[[i]]
        tr <- meta$Trait[[i]]
        ann_tbl <- get_annotation_data_list(ds, p_threshold = input$data_viz_gwas_threshold)[[tr]]
        hit_count <- if (is.null(ann_tbl)) 0 else nrow(ann_tbl)

        if (level == "Lipid") {
          return(tibble::tibble(Dataset = ds, Group = tr, Value = hit_count))
        }

        groups <- if (level == "Class") meta$ClassMembers[[i]] else meta$SubclassMembers[[i]]
        if (length(groups) == 0) {
          groups <- "Unclassified"
        }
        tibble::tibble(
          Dataset = ds,
          Group = groups,
          Value = hit_count / length(groups)
        )
      })
      joined <- bind_rows(rows)
      validate(need(nrow(joined) > 0, "No GWAS traits left after filtering."))
      joined
    }

    summary_df <- summary_df %>%
      mutate(Group = ifelse(is.na(Group) | Group == "", "Unclassified", Group)) %>%
      group_by(Dataset, Group) %>%
      summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop")

    top_n <- max(1, as.integer(input$data_viz_top_n))
    top_groups <- summary_df %>%
      group_by(Group) %>%
      summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(Total)) %>%
      head(top_n) %>%
      pull(Group)

    summary_df <- summary_df %>%
      mutate(Group = ifelse(Group %in% top_groups, Group, "Other")) %>%
      group_by(Dataset, Group) %>%
      summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop") %>%
      group_by(Dataset) %>%
      mutate(Percent = 100 * Value / sum(Value, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(
        DatasetLabel = ifelse(Dataset == "lowinput", "Low Input", "Control")
      )

    if (input$data_viz_mode == "single") {
      ds_label <- ifelse(datasets[[1]] == "lowinput", "Low Input", "Control")
      summary_df <- summary_df %>% mutate(DatasetLabel = ds_label)
    } else {
      summary_df$DatasetLabel <- factor(summary_df$DatasetLabel, levels = c("Low Input", "Control"))
    }

    group_order <- summary_df %>%
      group_by(Group) %>%
      summarise(Total = sum(Percent, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(Total)) %>%
      pull(Group)

    summary_df$Group <- factor(summary_df$Group, levels = group_order)
    palette <- grDevices::hcl.colors(length(group_order), palette = "Set 3")
    names(palette) <- group_order

    legend_labels <- if (input$data_viz_mode == "single") {
      lbl_tbl <- summary_df %>%
        group_by(Group) %>%
        summarise(PercentTotal = sum(Percent, na.rm = TRUE), .groups = "drop")
      setNames(
        sprintf("%s (%.1f%%)", as.character(lbl_tbl$Group), lbl_tbl$PercentTotal),
        as.character(lbl_tbl$Group)
      )
    } else {
      lbl_tbl <- summary_df %>%
        group_by(Group) %>%
        summarise(
          LowInput = sum(Percent[DatasetLabel == "Low Input"], na.rm = TRUE),
          Control = sum(Percent[DatasetLabel == "Control"], na.rm = TRUE),
          .groups = "drop"
        )
      setNames(
        sprintf("%s (L %.1f%%, C %.1f%%)", as.character(lbl_tbl$Group), lbl_tbl$LowInput, lbl_tbl$Control),
        as.character(lbl_tbl$Group)
      )
    }

    x_label <- if (input$data_viz_domain == "abundance") {
      "% Composition (mean across samples)"
    } else {
      "% Composition (GWAS hit counts)"
    }

    ggplot(summary_df, aes(x = Percent, y = DatasetLabel, fill = Group)) +
      geom_col(width = 0.6, color = "white", size = 0.2) +
      geom_text(
        data = summary_df %>% filter(Percent >= 6),
        aes(label = sprintf("%.1f%%", Percent)),
        position = position_stack(vjust = 0.5),
        size = 3,
        color = "black"
      ) +
      scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%")) +
      scale_fill_manual(values = palette, breaks = group_order, labels = legend_labels[group_order]) +
      labs(
        x = x_label,
        y = NULL,
        fill = level,
        title = paste("Composition by", level)
      ) +
      plot_theme +
      theme(
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank()
      )
  })

  output$data_visualize_summary <- renderUI({
    req(input$data_viz_domain, input$data_viz_mode, input$data_viz_level)
    datasets <- data_viz_selected_datasets()
    level_effective <- effective_data_viz_level()
    ds_label <- if (length(datasets) == 1) {
      ifelse(datasets[[1]] == "lowinput", "Low Input", "Control")
    } else {
      "Control + Low Input"
    }
    domain_label <- ifelse(input$data_viz_domain == "gwas", "GWAS traits", "Lipid abundance")
    tags$p(
      style = "margin-top: 10px; color: #555;",
      tags$strong("Mode: "), domain_label, "  |  ",
      tags$strong("Dataset: "), ds_label, "  |  ",
      tags$strong("Level: "), level_effective,
      if (input$data_viz_domain == "gwas") {
        tags$span("  |  ", tags$strong("Threshold: "), paste0("-log10(p) >= ", input$data_viz_gwas_threshold))
      },
      if (level_effective != input$data_viz_level) {
        tags$span("  |  Auto detail: showing subclass composition for the selected class.")
      }
    )
  })
  
  
  ##======================================================================
  ##  TAB 2: Boxplots with Sub-tabs (Individual, Class, Ratio Individual, Ratio Class)
  ##======================================================================
  
  # UI for lipid class dropdown (with "All")
  output$class_filter_plot_ui <- renderUI({
    classes <- sort(unique(na.omit(lipid_class_info$Class)))
    selectInput(
      inputId = "class_filter_plot",
      label = "Filter by Lipid Class:",
      choices = c("All", classes),
      selected = "All"
    )
  })
  
  # UI for lipid compound dropdown, includes "All (Sum)"
  output$compound_ui <- renderUI({
    req(input$plot_dataset)
    df <- dataset_list[[input$plot_dataset]]
    all_lipids <- colnames(df)[-1]
    
    if (!is.null(input$class_filter_plot) && input$class_filter_plot != "All") {
      selected_lipids <- lipid_class_info$Lipids[lipid_class_info$Class == input$class_filter_plot]
      available_lipids <- intersect(selected_lipids, all_lipids)
      choices <- c("All (Sum)" = "ALL_LIPIDS_SUM", available_lipids)
    } else {
      available_lipids <- all_lipids
      choices <- available_lipids
    }
    
    selectInput(
      inputId = "compound_choice",
      label = "Choose a Lipid Compound:",
      choices = choices,
      selected = available_lipids[1]
    )
  })
  
  # Class filter
  output$class_filter_plot_ui <- renderUI({
    classes <- sort(unique(na.omit(lipid_class_info$Class)))
    selectInput("class_filter_plot", "Filter by Lipid Class:", choices = c("All", classes), selected = "All")
  })
  
  # Subclass filter (conditional)
  output$subclass_filter_ui <- renderUI({
    req(input$class_filter_plot)
    if (input$class_filter_plot == "All") return(NULL)
    
    subclasses <- lipid_class_info %>%
      filter(Class == input$class_filter_plot, !is.na(SubClass)) %>%
      distinct(SubClass) %>%
      pull(SubClass)
    
    if (length(subclasses) == 0) return(NULL)
    
    selectInput("subclass_filter", "Filter by Subclass:", choices = c("All", subclasses), selected = "All")
  })
  
  # Sub-subclass filter (conditional)
  output$subsubclass_filter_ui <- renderUI({
    req(input$class_filter_plot, input$subclass_filter)
    if (input$class_filter_plot == "All" || input$subclass_filter == "All") return(NULL)
    if (is.na(subsubclass_col)) return(NULL)
    
    subsubclasses <- lipid_class_info %>%
      filter(
        Class == input$class_filter_plot,
        SubClass == input$subclass_filter,
        !is.na(.data[[subsubclass_col]])
      ) %>%
      distinct(.data[[subsubclass_col]]) %>%
      pull(.data[[subsubclass_col]])
    
    if (length(subsubclasses) == 0) return(NULL)
    
    selectInput("subsubclass_filter", "Filter by Sub-subclass:", choices = c("All", subsubclasses), selected = "All")
  })
  
  # Compound filter
  output$compound_ui <- renderUI({
    req(input$plot_dataset)
    df <- dataset_list[[input$plot_dataset]]
    all_lipids <- colnames(df)[-1]
    
    filtered_lipids <- lipid_class_info
    
    if (!is.null(input$class_filter_plot) && input$class_filter_plot != "All") {
      filtered_lipids <- filtered_lipids %>% filter(Class == input$class_filter_plot)
    }
    
    if (!is.null(input$subclass_filter) && input$subclass_filter != "All") {
      filtered_lipids <- filtered_lipids %>% filter(SubClass == input$subclass_filter)
    }
    
    if (!is.na(subsubclass_col) && !is.null(input$subsubclass_filter) && input$subsubclass_filter != "All") {
      filtered_lipids <- filtered_lipids %>% filter(.data[[subsubclass_col]] == input$subsubclass_filter)
    }
    
    valid_lipids <- intersect(filtered_lipids$Lipids, all_lipids)
    if (length(valid_lipids) == 0) {
      return(selectInput("compound_choice", "Choose a Lipid Compound:", choices = character(0)))
    }
    choices <- c("All (Sum)" = "ALL_LIPIDS_SUM", valid_lipids)
    
    selectInput("compound_choice", "Choose a Lipid Compound:", choices = choices, selected = valid_lipids[1])
  })
  
  # # Single dataset plot
  output$boxplot_single <- renderPlot({
    req(input$plot_dataset, input$compound_choice, input$data_type_boxplot)
    df <- get_dataset(input$plot_dataset, input$data_type_boxplot)

    if (input$compound_choice == "ALL_LIPIDS_SUM") {
      filtered_lipids <- lipid_class_info

      if (!is.null(input$class_filter_plot) && input$class_filter_plot != "All") {
        filtered_lipids <- filtered_lipids %>% filter(Class == input$class_filter_plot)
      }

      if (!is.null(input$subclass_filter) && input$subclass_filter != "All") {
        filtered_lipids <- filtered_lipids %>% filter(SubClass == input$subclass_filter)
      }

      if (!is.na(subsubclass_col) && !is.null(input$subsubclass_filter) && input$subsubclass_filter != "All") {
        filtered_lipids <- filtered_lipids %>% filter(.data[[subsubclass_col]] == input$subsubclass_filter)
      }

      valid_lipids <- intersect(filtered_lipids$Lipids, colnames(df))
      validate(need(length(valid_lipids) > 0, "No lipids available for the selected class/subclass filters."))

      df_plot <- df %>%
        select(Sample = 1, all_of(valid_lipids)) %>%
        mutate(Value = rowSums(across(all_of(valid_lipids)), na.rm = TRUE)) %>%
        select(Sample, Value) %>%
        mutate(Condition = input$plot_dataset)

      scope_label <- "selected lipids"
      if (!is.na(subsubclass_col) && !is.null(input$subsubclass_filter) && input$subsubclass_filter != "All") {
        scope_label <- paste(input$subsubclass_filter, "sub-subclass")
      } else if (!is.null(input$subclass_filter) && input$subclass_filter != "All") {
        scope_label <- paste(input$subclass_filter, "subclass")
      } else if (!is.null(input$class_filter_plot) && input$class_filter_plot != "All") {
        scope_label <- paste(input$class_filter_plot, "class")
      }

      plot_title <- paste("Summed abundance of", scope_label)

    } else {
      df_plot <- df %>%
        select(Sample = 1, Value = all_of(input$compound_choice)) %>%
        mutate(Condition = input$plot_dataset)

      plot_title <- paste("Distribution of", input$compound_choice, "in", input$plot_dataset)
    }
    
    validate(
      need(nrow(df_plot) > 1, "Not enough data points to draw this plot."),
      need(any(is.finite(df_plot$Value)), "No finite values available for this selection.")
    )

    # Dynamic x-axis label based on data type
    x_label <- if (input$data_type_boxplot == "raw") "Raw Signal Intensity" else data_type_label(input$data_type_boxplot)

    p <- ggplot(df_plot, aes(x = Value)) +
      geom_density(
        fill = ifelse(input$plot_dataset == "control", "#FF7F0E", "#2CA02C"),
        alpha = 0.6,
        color = NA
      ) +
      geom_vline(aes(xintercept = mean(Value, na.rm = TRUE), color = "Mean"),
                 linetype = "dashed", size = 1, show.legend = TRUE) +
      geom_vline(aes(xintercept = median(Value, na.rm = TRUE), color = "Median"),
                 linetype = "solid", size = 1, show.legend = TRUE) +
      scale_fill_manual(values = c("control" = "#1f77b4", "lowinput" = "#ff7f0e")) +
      scale_color_manual(values = c("Mean" = "blue", "Median" = "red")) +
      labs(
        title = plot_title,
        x = x_label,
        y = "Density",
        fill = NULL,
        color = NULL
      ) +
      plot_theme

    plot_store$boxplot_single <- p
    p
  })

  # Combined dataset plot
  output$boxplot_both <- renderPlot({
    req(input$compound_choice, input$plot_dataset, input$class_filter_plot, input$data_type_boxplot)

    other_dataset <- setdiff(names(dataset_list), input$plot_dataset)
    df_main <- get_dataset(input$plot_dataset, input$data_type_boxplot)
    df_other <- get_dataset(other_dataset, input$data_type_boxplot)

    if (input$compound_choice == "ALL_LIPIDS_SUM") {
      filtered_lipids <- lipid_class_info

      if (!is.null(input$class_filter_plot) && input$class_filter_plot != "All") {
        filtered_lipids <- filtered_lipids %>% filter(Class == input$class_filter_plot)
      }

      if (!is.null(input$subclass_filter) && input$subclass_filter != "All") {
        filtered_lipids <- filtered_lipids %>% filter(SubClass == input$subclass_filter)
      }

      if (!is.na(subsubclass_col) && !is.null(input$subsubclass_filter) && input$subsubclass_filter != "All") {
        filtered_lipids <- filtered_lipids %>% filter(.data[[subsubclass_col]] == input$subsubclass_filter)
      }

      valid_main <- intersect(filtered_lipids$Lipids, colnames(df_main))
      valid_other <- intersect(filtered_lipids$Lipids, colnames(df_other))

      if (length(valid_main) == 0 || length(valid_other) == 0) return(NULL)

      df_main_sum <- df_main %>%
        dplyr::select(Sample = 1, all_of(valid_main)) %>%
        dplyr::mutate(Value = rowSums(across(all_of(valid_main)), na.rm = TRUE),
               Condition = input$plot_dataset) %>%
        dplyr::select(Sample, Value, Condition)

      df_other_sum <- df_other %>%
        dplyr::select(Sample = 1, all_of(valid_other)) %>%
        dplyr::mutate(Value = rowSums(across(all_of(valid_other)), na.rm = TRUE),
               Condition = other_dataset) %>%
        dplyr::select(Sample, Value, Condition)

      df_combined <- bind_rows(df_main_sum, df_other_sum)

      scope_label <- "selected lipids"
      if (!is.na(subsubclass_col) && !is.null(input$subsubclass_filter) && input$subsubclass_filter != "All") {
        scope_label <- paste(input$subsubclass_filter, "sub-subclass")
      } else if (!is.null(input$subclass_filter) && input$subclass_filter != "All") {
        scope_label <- paste(input$subclass_filter, "subclass")
      } else if (!is.null(input$class_filter_plot) && input$class_filter_plot != "All") {
        scope_label <- paste(input$class_filter_plot, "class")
      }

      plot_title <- paste("Summed abundance of", scope_label, "across conditions")

    } else {
      compound_col <- input$compound_choice

      if (!(compound_col %in% colnames(df_main)) || !(compound_col %in% colnames(df_other))) {
        return(NULL)
      }

      df_main_plot <- df_main %>%
        select(Sample = 1, Value = all_of(compound_col)) %>%
        mutate(Condition = input$plot_dataset)

      df_other_plot <- df_other %>%
        select(Sample = 1, Value = all_of(compound_col)) %>%
        mutate(Condition = other_dataset)

      df_combined <- bind_rows(df_main_plot, df_other_plot)

      plot_title <- paste("Comparison of", compound_col, "across conditions")
    }
    
    validate(
      need(nrow(df_combined) > 1, "Not enough data points to compare conditions."),
      need(any(is.finite(df_combined$Value)), "No finite values available for this selection.")
    )

    # Dynamic y-axis label based on data type
    y_label <- if (input$data_type_boxplot == "raw") "Raw Signal Intensity" else data_type_label(input$data_type_boxplot)

    # Statistical test (Wilcoxon)
    ctrl_vals <- df_combined$Value[df_combined$Condition == "control"]
    low_vals <- df_combined$Value[df_combined$Condition == "lowinput"]

    pval <- tryCatch({
      wilcox.test(ctrl_vals, low_vals)$p.value
    }, error = function(e) NA)

    # Format p-value for display
    if (!is.na(pval)) {
      if (pval < 0.001) {
        pval_text <- "p < 0.001"
      } else {
        pval_text <- paste0("p = ", signif(pval, 3))
      }
      significance <- ifelse(pval < 0.05, "*", "ns")
    } else {
      pval_text <- ""
      significance <- ""
    }

    # Calculate y position for annotation
    y_max <- max(df_combined$Value, na.rm = TRUE)

    p <- ggplot(df_combined, aes(x = Condition, y = Value, fill = Condition)) +
      geom_violin(alpha = 0.5, width = 0.8) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      scale_fill_manual(values = c("control" = "#FF7F0E", "lowinput" = "#2CA02C")) +
      labs(
        title = plot_title,
        subtitle = if (pval_text != "") paste0("Wilcoxon test: ", pval_text, " ", significance) else NULL,
        x = "Condition",
        y = y_label
      ) +
      plot_theme +
      theme(legend.position = "none",
            plot.subtitle = element_text(size = 12, color = ifelse(!is.na(pval) && pval < 0.05, "#27ae60", "#666")))

    # Add significance bracket if significant
    if (!is.na(pval) && pval < 0.05) {
      p <- p +
        annotate("segment", x = 1, xend = 2, y = y_max * 1.05, yend = y_max * 1.05) +
        annotate("text", x = 1.5, y = y_max * 1.1, label = significance, size = 6)
    }

    plot_store$boxplot_both <- p
    p
  })
  
  #=======================
  # Ratio Tabs (UI)
  #=======================
  
  output$ratio_ui_1 <- renderUI({
    req(input$plot_dataset)
    df <- dataset_list[[input$plot_dataset]]
    all_lipids <- colnames(df)[-1]
    sum_lipids <- paste0("SUM_", sort(unique(na.omit(lipid_class_info$Class))))
    selectInput("ratio_lipid_1", "Numerator Lipid:", choices = c(sum_lipids, all_lipids))
  })
  
  output$ratio_ui_2 <- renderUI({
    req(input$plot_dataset)
    df <- dataset_list[[input$plot_dataset]]
    all_lipids <- colnames(df)[-1]
    sum_lipids <- paste0("SUM_", sort(unique(na.omit(lipid_class_info$Class))))
    selectInput("ratio_lipid_2", "Denominator Lipid:", choices = c(sum_lipids, all_lipids))
  })
  
  #=======================
  # Shared Function for Boxplot Rendering
  #=======================
  
  plot_lipid_boxplot <- function(df, values, title, condition) {
    df_plot <- data.frame(Sample = df[[1]], Value = values, Condition = condition)

    ggplot(df_plot, aes(x = Value)) +
      geom_density(fill = ifelse(condition == "control", "#FF7F0E", "#2CA02C"), alpha = 0.6) +
      geom_vline(aes(xintercept = mean(Value, na.rm = TRUE), color = "Mean"),
                 linetype = "dashed", linewidth = 1, show.legend = TRUE) +
      geom_vline(aes(xintercept = median(Value, na.rm = TRUE), color = "Median"),
                 linetype = "solid", linewidth = 1, show.legend = TRUE) +
      scale_color_manual(values = c("Mean" = "blue", "Median" = "red")) +
      labs(
        title = title,
        x = "log10(Ratio)",
        y = "Density",
        color = ""
      ) +
      plot_theme
  }
  
  get_ratio_values <- function(df, lipid_input) {
    if (startsWith(lipid_input, "SUM_")) {
      class_name <- sub("^SUM_", "", lipid_input)
      lipids <- lipid_class_info$Lipids[lipid_class_info$Class == class_name]
      valid <- intersect(lipids, colnames(df))
      if (length(valid) == 0) return(rep(NA, nrow(df)))  # guard clause
      rowSums(df[, valid, drop = FALSE], na.rm = TRUE)
    } else {
      df[[lipid_input]]
    }
  }
  
  #=======================
  # Ratio Individual Dataset Plot
  #=======================
  
  output$boxplot_ratio <- renderPlot({
    req(input$plot_dataset, input$ratio_lipid_1, input$ratio_lipid_2)
    df <- dataset_list[[input$plot_dataset]]

    numerator <- get_ratio_values(df, input$ratio_lipid_1)
    denominator <- get_ratio_values(df, input$ratio_lipid_2)

    if (any(is.na(numerator)) || any(is.na(denominator))) return(NULL)

    ratio_vals <- log10(numerator / denominator)

    p <- plot_lipid_boxplot(df, ratio_vals,
                       paste("log10(Ratio):", input$ratio_lipid_1, "/", input$ratio_lipid_2),
                       input$plot_dataset)
    plot_store$boxplot_ratio <- p
    p
  })
  
  #=======================
  # Ratio Combined Plot (Control vs LowInput)
  #=======================
  
  output$boxplot_ratio_combined <- renderPlot({
    req(input$ratio_lipid_1, input$ratio_lipid_2)
    datasets <- names(dataset_list)
    
    df_list <- lapply(datasets, function(ds) {
      df <- dataset_list[[ds]]
      numerator <- get_ratio_values(df, input$ratio_lipid_1)
      denominator <- get_ratio_values(df, input$ratio_lipid_2)
      
      if (any(is.na(numerator)) || any(is.na(denominator))) return(NULL)
      
      ratio_vals <- log10(numerator / denominator)
      data.frame(Sample = df[[1]], Value = ratio_vals, Condition = ds)
    })
    
    df_combined <- do.call(rbind, df_list[!sapply(df_list, is.null)])
    
    # ✅ Make sure there are exactly two levels
    if (length(unique(df_combined$Condition)) < 2) return(NULL)
    
    # Perform t-test
    ttest <- t.test(Value ~ Condition, data = df_combined)
    pval <- ttest$p.value
    
    # Format significance
    sig_label <- if (pval < 0.001) {
      "***"
    } else if (pval < 0.01) {
      "**"
    } else if (pval < 0.05) {
      "*"
    } else {
      "ns"
    }
    
    # 🔥 Trim for plot only
    lower <- quantile(df_combined$Value, 0.025, na.rm = TRUE)
    upper <- quantile(df_combined$Value, 0.975, na.rm = TRUE)
    df_trimmed <- df_combined %>% filter(Value >= lower & Value <= upper)
    
    p <- ggplot(df_trimmed, aes(x = Condition, y = Value, fill = Condition)) +
      geom_violin(alpha = 0.5, width = 0.8) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      geom_text(
        aes(x = 1.5, y = max(df_trimmed$Value, na.rm = TRUE) * 1.05, label = sig_label),
        inherit.aes = FALSE,
        size = 6
      ) +
      scale_fill_manual(values = c("control" = "#FF7F0E", "lowinput" = "#2CA02C")) +
      labs(
        title = paste("log10(Ratio of", input$ratio_lipid_1, "/", input$ratio_lipid_2, ") across conditions"),
        y = "log10(Ratio)",
        x = "Condition"
      ) +
      plot_theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )

    plot_store$boxplot_ratio_combined <- p
    p
  })
  
  
  ##======================================================================
  ##  TAB 3: PCA
  ##======================================================================
 
  output$pca_selected <- renderPrint({
    eventdata <- event_data("plotly_selected")
    if (is.null(eventdata)) {
      cat("Select points on the plot to see their labels here.")
    } else {
      unique(eventdata$text)
    }
  })
  
  output$pca_plot <- renderPlotly({  # plotly for interaction
    req(input$pca_dataset, input$data_type_pca)
    pca_data_mode <- if (input$data_type_pca == "clr") "tic" else input$data_type_pca
    df <- get_dataset(input$pca_dataset, pca_data_mode)
    pca_bundle <- compute_pca_bundle(
      df = df,
      mode = input$pca_mode,
      data_type = input$data_type_pca,
      class_filter = if (input$pca_mode == "Lipids") input$pca_class_filter else "All",
      center = input$pca_center,
      scale = input$pca_scale
    )
    validate(need(!is.null(pca_bundle), "PCA not available for this selection."))

    pca_data <- pca_bundle$scores
    xlab <- sprintf("PC1 (%.1f%%)", 100 * pca_bundle$pve[1])
    ylab <- sprintf("PC2 (%.1f%%)", 100 * pca_bundle$pve[2])

    if (input$pca_mode == "Samples") {
      plot_ly(
        data = pca_data,
        x = ~PC1,
        y = ~PC2,
        text = ~ID,
        type = "scatter",
        mode = "markers",
        marker = list(size = 7)
      ) %>%
        layout(
          title = paste("PCA - Samples (", input$pca_dataset, ",", pca_bundle$data_label, ")"),
          xaxis = list(title = xlab),
          yaxis = list(title = ylab),
          dragmode = "select"
        )
    } else {
      p <- build_lipid_pca_plot(
        pca_data = pca_data,
        title_text = paste("PCA - Lipids (", input$pca_dataset, ",", pca_bundle$data_label, ")"),
        xlab = xlab,
        ylab = ylab
      )

      ggplotly(p, tooltip = "text")
    }
  })
  
  
  ##======================================================================
  ##  TAB 4: t-SNE / UMAP
  ##======================================================================
  
  # output$dr_selected <- renderPrint({
  #   eventdata <- event_data("plotly_selected")
  #   if (is.null(eventdata)) {
  #     cat("Select points on the plot to see labels here.")
  #   } else {
  #     unique(eventdata$text)
  #   }
  # })
  # 
  # output$dr_plot <- renderPlotly({
  #   df <- dataset_list[[input$dr_dataset]]
  #   method <- input$dr_method
  #   
  #   # Handle sample vs lipid mode
  #   if (input$dr_mode == "Samples") {
  #     numeric_df <- df[, -1]
  #     numeric_df <- as.data.frame(sapply(numeric_df, as.numeric))
  #     row_labels <- df[[1]]
  #   } else {
  #     numeric_df <- df[, -1]
  #     numeric_df <- t(numeric_df)
  #     colnames(numeric_df) <- df[[1]]
  #     numeric_df <- as.data.frame(numeric_df)
  #     row_labels <- rownames(numeric_df)
  #   }
  #   
  #   # Dimensionality reduction
  #   if (method == "t-SNE") {
  #     library(Rtsne)
  #     dr_out <- Rtsne(numeric_df, perplexity = input$perplexity, check_duplicates = FALSE)
  #     dim_df <- data.frame(Dim1 = dr_out$Y[,1], Dim2 = dr_out$Y[,2], Label = row_labels)
  #   } else {
  #     library(umap)
  #     dr_out <- umap(numeric_df)
  #     dim_df <- data.frame(Dim1 = dr_out$layout[,1], Dim2 = dr_out$layout[,2], Label = row_labels)
  #   }
  #   
  #   plot_ly(
  #     data = dim_df,
  #     x = ~Dim1,
  #     y = ~Dim2,
  #     type = "scatter",
  #     mode = "markers",
  #     text = ~Label,
  #     marker = list(size = 6)
  #   ) %>%
  #     layout(
  #       title = paste(method, "-", input$dr_mode, "(", input$dr_dataset, ")"),
  #       xaxis = list(title = "Dim 1"),
  #       yaxis = list(title = "Dim 2"),
  #       dragmode = "select"
  #     )
  # })
  
  output$dr_selected <- renderPrint({
    eventdata <- event_data("plotly_selected")
    if (is.null(eventdata)) {
      cat("Select points on the plot to see labels here.")
    } else {
      unique(eventdata$text)
    }
  })
  
  output$dr_plot <- renderPlotly({
    req(input$dr_dataset, input$data_type_dr)
    df <- get_dataset(input$dr_dataset, input$data_type_dr)
    method <- input$dr_method

    if (input$dr_mode == "Samples") {
      numeric_df <- df[, -1]
      numeric_df <- as.data.frame(sapply(numeric_df, as.numeric))
      row_labels <- df[[1]]

      # Dimensionality reduction
      if (method == "t-SNE") {
        library(Rtsne)
        dr_out <- Rtsne(numeric_df, perplexity = input$perplexity, check_duplicates = FALSE)
        dim_df <- data.frame(Dim1 = dr_out$Y[,1], Dim2 = dr_out$Y[,2], Label = row_labels)
      } else {
        library(umap)
        dr_out <- umap(numeric_df)
        dim_df <- data.frame(Dim1 = dr_out$layout[,1], Dim2 = dr_out$layout[,2], Label = row_labels)
      }

      # Data type label for title
      data_label <- data_type_label(input$data_type_dr)

      # Simple plotly for Samples mode
      plot_ly(
        data = dim_df,
        x = ~Dim1,
        y = ~Dim2,
        type = "scatter",
        mode = "markers",
        text = ~Label,
        marker = list(size = 6)
      ) %>%
        layout(
          title = paste(method, "- Samples (", input$dr_dataset, ",", data_label, ")"),
          xaxis = list(title = "Dim 1"),
          yaxis = list(title = "Dim 2"),
          dragmode = "select"
        )

    } else {
      # Lipids mode - with ellipses
      numeric_df <- df[, -1]
      numeric_df <- t(numeric_df)
      colnames(numeric_df) <- df[[1]]
      numeric_df <- as.data.frame(numeric_df)
      row_labels <- rownames(numeric_df)

      # Extract class from lipid name using get_lipid_class (like paper script)
      class_col <- get_lipid_class(row_labels)

      # Dimensionality reduction
      if (method == "t-SNE") {
        dr_out <- Rtsne(numeric_df, perplexity = input$perplexity, check_duplicates = FALSE)
        dim_df <- data.frame(Dim1 = dr_out$Y[,1], Dim2 = dr_out$Y[,2], Label = row_labels, Class = class_col)
      } else {
        dr_out <- umap(numeric_df)
        dim_df <- data.frame(Dim1 = dr_out$layout[,1], Dim2 = dr_out$layout[,2], Label = row_labels, Class = class_col)
      }

      # Add Group information
      dim_df$Group <- class_to_group(dim_df$Class)

      # Filter out "Other" class lipids for cleaner visualization
      dim_df <- dim_df[dim_df$Class != "Other", ]
      # Also filter out NA groups
      dim_df <- dim_df[!is.na(dim_df$Group), ]

      # Data type label for title
      data_label <- data_type_label(input$data_type_dr)

      # ggplot with ellipses (exactly like paper script)
      # Points colored by Class, ellipses by Group
      p <- ggplot(dim_df, aes(x = Dim1, y = Dim2)) +
        stat_ellipse(
          aes(color = Group),
          level = 0.95,
          linewidth = 1.5,
          type = "norm",
          show.legend = FALSE
        ) +
        scale_color_manual(values = group_colors) +
        ggnewscale::new_scale_color() +
        geom_point(aes(color = Class, text = Label), size = 3, alpha = 0.9) +
        scale_color_manual(values = class_colors, name = "Class") +
        labs(
          title = paste(method, "- Lipids (", input$dr_dataset, ",", data_label, ")"),
          x = "Dim 1",
          y = "Dim 2"
        ) +
        plot_theme +
        theme(legend.position = "right")

      ggplotly(p, tooltip = "text")
    }
  })
  
  ##======================================================================
  ##  TAB 5: Volcano Plot
  ##======================================================================
  output$volcano_plot <- renderPlot({
    # We'll do a simple approach:
    # 1) Find columns in both 'control' and 'lowinput'
    # 2) For each of those columns, run a t-test
    # 3) log2FC = mean(lowinput) - mean(control) in log2
    #    pval from t-test
    
    common_cols <- intersect(colnames(control), colnames(lowinput))
    # remove the first column (Sample ID) from analysis
    common_cols <- setdiff(common_cols, colnames(control)[1])
    
    # if no common lipid columns, just bail
    validate(
      need(length(common_cols) > 0, "No overlapping lipid columns found!")
    )
    
    # Build a results data frame
    results <- data.frame(Lipid = common_cols, log2FC = NA, pval = NA)
    
    for (i in seq_along(common_cols)) {
      col_i <- common_cols[i]
      # Extract numeric vectors
      v_control  <- as.numeric(control[[col_i]])
      v_lowinput <- as.numeric(lowinput[[col_i]])
      
      # T-test
      # For fold change: typically (LowInput / Control). We'll do log2 of means
      t_out <- t.test(v_lowinput, v_control)
      mean_control  <- mean(v_control,  na.rm = TRUE)
      mean_lowinput <- mean(v_lowinput, na.rm = TRUE)
      # Avoid zero or negative means in log (in real data handle carefully)
      # We'll do log2FC = log2( mean_lowinput / mean_control )
      fc       = mean_lowinput / mean_control
      log2fc   = log2(fc)
      
      results$log2FC[i] <- log2fc
      results$pval[i]   <- t_out$p.value
    }
    
    # -log10 p-value
    results$negLog10p <- -log10(results$pval)
    
    # Basic volcano
    ggplot(results, aes(x = log2FC, y = negLog10p)) +
      geom_point() +
      labs(
        title = "Volcano Plot (LowInput vs Control)",
        x = "log2 Fold Change",
        y = "-log10 p-value"
      ) +
      plot_theme
  })
  
  ##======================================================================
  ##  TAB 6: Heatmap
  ##======================================================================

  # UI for individual lipid selection in heatmap
  output$heatmap_lipid_selector <- renderUI({
    req(input$heatmap_dataset, input$data_type_heatmap)

    if (input$heatmap_dataset == "Both") {
      df_control <- get_dataset("control", input$data_type_heatmap)
      df_lowinput <- get_dataset("lowinput", input$data_type_heatmap)
      all_lipids <- sort(union(colnames(df_control)[-1], colnames(df_lowinput)[-1]))
    } else {
      df <- get_dataset(input$heatmap_dataset, input$data_type_heatmap)
      all_lipids <- colnames(df)[-1]
    }

    selectizeInput("heatmap_selected_lipids", "Select Lipids:",
                   choices = all_lipids,
                   selected = all_lipids[1:min(10, length(all_lipids))],
                   multiple = TRUE,
                   options = list(maxItems = 50, placeholder = "Select lipids..."))
  })

  output$heatmap_plot <- renderPlotly({
    req(input$heatmap_dataset, input$data_type_heatmap, input$heatmap_lipid_selection)

    # Handle "Both" case - combine control and lowinput with unique lipids
    if (input$heatmap_dataset == "Both") {
      df_control <- get_dataset("control", input$data_type_heatmap)
      df_lowinput <- get_dataset("lowinput", input$data_type_heatmap)

      # Get all lipids (union)
      all_lipids <- union(colnames(df_control)[-1], colnames(df_lowinput)[-1])

      # Create matrices with NA for missing lipids
      ctrl_mat <- df_control[, -1]
      ctrl_mat <- as.data.frame(sapply(ctrl_mat, as.numeric))
      rownames(ctrl_mat) <- paste0("ctrl_", df_control[[1]])

      low_mat <- df_lowinput[, -1]
      low_mat <- as.data.frame(sapply(low_mat, as.numeric))
      rownames(low_mat) <- paste0("low_", df_lowinput[[1]])

      # Add missing columns with NA
      missing_in_ctrl <- setdiff(all_lipids, colnames(ctrl_mat))
      missing_in_low <- setdiff(all_lipids, colnames(low_mat))

      for (lip in missing_in_ctrl) {
        ctrl_mat[[lip]] <- NA_real_
      }
      for (lip in missing_in_low) {
        low_mat[[lip]] <- NA_real_
      }

      # Reorder columns to match
      ctrl_mat <- ctrl_mat[, all_lipids, drop = FALSE]
      low_mat <- low_mat[, all_lipids, drop = FALSE]

      # Combine both matrices
      mat <- rbind(ctrl_mat, low_mat)

    } else {
      df <- get_dataset(input$heatmap_dataset, input$data_type_heatmap)

      # Strip first column and ensure numeric
      mat <- df[, -1]
      mat <- as.data.frame(sapply(mat, as.numeric))
      rownames(mat) <- df[[1]]
    }

    # Remove completely NA rows/columns
    mat <- mat[rowSums(is.na(mat)) != ncol(mat), ]
    mat <- mat[, colSums(is.na(mat)) != nrow(mat)]

    # Direction of comparison: Samples vs Lipids
    axis_mode <- ifelse(input$heatmap_mode == "Samples", 1, 2)

    # Select lipids based on selection mode
    lipid_selection <- input$heatmap_lipid_selection

    if (lipid_selection == "top_changed") {
      # Mean absolute deviation from median (handle NA)
      diff_metric <- apply(mat, axis_mode, function(x) {
        med <- median(x, na.rm = TRUE)
        mean(abs(x - med), na.rm = TRUE)
      })
      top_n <- input$top_n_heatmap
      selected <- names(sort(diff_metric, decreasing = TRUE))[1:min(top_n, length(diff_metric))]

      # Subset matrix
      if (axis_mode == 1) {
        mat_top <- mat[selected, , drop = FALSE]
      } else {
        mat_top <- mat[, selected, drop = FALSE]
      }

    } else if (lipid_selection == "by_class") {
      # Filter by lipid class
      req(input$heatmap_class_filter)
      lipid_classes <- sapply(colnames(mat), get_lipid_class)
      selected_lipids <- names(lipid_classes)[lipid_classes == input$heatmap_class_filter]
      selected_lipids <- intersect(selected_lipids, colnames(mat))

      validate(need(length(selected_lipids) > 1, paste("Not enough lipids found for class:", input$heatmap_class_filter)))

      mat_top <- mat[, selected_lipids, drop = FALSE]

    } else if (lipid_selection == "individual") {
      # User-selected individual lipids
      req(input$heatmap_selected_lipids)
      selected_lipids <- intersect(input$heatmap_selected_lipids, colnames(mat))

      validate(need(length(selected_lipids) > 1, "Please select at least 2 lipids."))

      mat_top <- mat[, selected_lipids, drop = FALSE]
    }

    # Generate short names for columns (i.e., lipids)
    original_lipid_names <- colnames(mat_top)
    short_lipid_ids <- paste0("L", seq_along(original_lipid_names))
    lipid_name_map <- data.frame(ID = short_lipid_ids, Name = original_lipid_names)

    # Save mapping so you can render it later
    output$lipid_name_lookup <- DT::renderDataTable({
      lipid_name_map
    }, options = list(pageLength = 10, scrollX = TRUE))

    # Rename columns in the heatmap matrix
    colnames(mat_top) <- short_lipid_ids

    # Validate size
    validate(
      need(nrow(mat_top) > 1 && ncol(mat_top) > 1, "Not enough data for heatmap.")
    )

    # Scale the matrix (handling NA values)
    mat_scaled <- scale(mat_top)

    # Custom color palette with dark gray for NA values
    na_color <- "#2d2d2d"  # Dark gray/black for missing values

    # Get color palette based on user selection
    palette_choice <- input$heatmap_color_palette
    if (is.null(palette_choice)) palette_choice <- "RdYlBu"

    if (palette_choice == "viridis") {
      main_colors <- viridisLite::viridis(255)
    } else if (palette_choice == "plasma") {
      main_colors <- viridisLite::plasma(255)
    } else {
      main_colors <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(255)
    }

    # Plot
    heatmaply(
      mat_scaled,
      scale = "none",
      dendrogram = "both",
      showticklabels = c(TRUE, TRUE),
      colors = main_colors,
      na.value = na_color,
      xlab = ifelse(axis_mode == 1, "Lipids", "Lipids"),
      ylab = ifelse(axis_mode == 1, "Samples", "Samples"),
      labCol = NULL,
      labRow = rownames(mat_top),
      margins = c(40, 120, 40, 120)
    )
  })
  
  
  # UI for individual lipid selection in condition comparison heatmap
  output$heatmap_diff_lipid_selector <- renderUI({
    req(input$data_type_heatmap_diff)
    df_control <- get_dataset("control", input$data_type_heatmap_diff)
    df_lowinput <- get_dataset("lowinput", input$data_type_heatmap_diff)
    # Use ALL lipids (union) instead of just common ones
    all_lipids <- sort(union(colnames(df_control)[-1], colnames(df_lowinput)[-1]))
    selectizeInput("heatmap_diff_selected_lipids", "Select Lipids:",
                   choices = all_lipids, selected = head(all_lipids, 10),
                   multiple = TRUE, options = list(maxItems = 50))
  })

  output$heatmap_diff <- renderPlotly({
    req(input$data_type_heatmap_diff, input$heatmap_diff_lipid_selection)
    df_control <- get_dataset("control", input$data_type_heatmap_diff)
    df_lowinput <- get_dataset("lowinput", input$data_type_heatmap_diff)

    # Get numeric matrix
    ctrl_mat <- df_control[, -1]
    low_mat <- df_lowinput[, -1]
    rownames(ctrl_mat) <- df_control[[1]]
    rownames(low_mat) <- df_lowinput[[1]]
    ctrl_mat <- as.data.frame(sapply(ctrl_mat, as.numeric))
    low_mat <- as.data.frame(sapply(low_mat, as.numeric))

    # Get ALL lipids (union of both conditions) to include unique ones
    all_cols <- union(colnames(ctrl_mat), colnames(low_mat))
    common_cols <- intersect(colnames(ctrl_mat), colnames(low_mat))
    ctrl_only <- setdiff(colnames(ctrl_mat), colnames(low_mat))
    low_only <- setdiff(colnames(low_mat), colnames(ctrl_mat))

    # Z-score each condition by column (only for columns that exist)
    z_control <- scale(ctrl_mat)
    z_lowinput <- scale(low_mat)

    # Mean Z-score per lipid for each condition
    mean_z_control <- colMeans(z_control, na.rm = TRUE)
    mean_z_lowinput <- colMeans(z_lowinput, na.rm = TRUE)

    # Create full vectors with NA for missing lipids
    full_mean_z_control <- setNames(rep(NA_real_, length(all_cols)), all_cols)
    full_mean_z_lowinput <- setNames(rep(NA_real_, length(all_cols)), all_cols)
    full_mean_z_control[names(mean_z_control)] <- mean_z_control
    full_mean_z_lowinput[names(mean_z_lowinput)] <- mean_z_lowinput

    # Calculate difference (will be NA for unique lipids)
    z_diff <- full_mean_z_lowinput - full_mean_z_control

    # Get lipid selection mode
    lipid_selection <- input$heatmap_diff_lipid_selection

    if (lipid_selection == "top_changed") {
      # Top N most changed lipids (from common lipids only for ranking)
      top_n <- input$top_n_diff
      # For top changed, include both common lipids ranked by difference AND unique lipids
      common_z_diff <- z_diff[common_cols]
      ranked_common <- names(sort(abs(common_z_diff), decreasing = TRUE))
      # Take top N from common, then add unique lipids
      top_common <- head(ranked_common, top_n)
      # Include some unique lipids if they exist
      top_lipids <- c(top_common, ctrl_only, low_only)
      top_lipids <- head(top_lipids, top_n + length(ctrl_only) + length(low_only))
    } else if (lipid_selection == "by_class") {
      # Filter by lipid class - include ALL lipids of that class
      req(input$heatmap_diff_class_filter)
      lipid_classes <- sapply(all_cols, get_lipid_class)
      top_lipids <- names(lipid_classes)[lipid_classes == input$heatmap_diff_class_filter]
      validate(need(length(top_lipids) > 1, paste("Not enough lipids found for class:", input$heatmap_diff_class_filter)))
    } else if (lipid_selection == "individual") {
      # Individual lipid selection
      req(input$heatmap_diff_selected_lipids)
      top_lipids <- intersect(input$heatmap_diff_selected_lipids, all_cols)
      validate(need(length(top_lipids) > 1, "Please select at least 2 lipids."))
    }

    # Subset and make matrix
    mat_comp <- rbind(full_mean_z_control[top_lipids], full_mean_z_lowinput[top_lipids])
    rownames(mat_comp) <- c("control", "lowinput")

    # Generate short names for columns (i.e., lipids) - same as Individual Samples tab
    original_lipid_names <- colnames(mat_comp)
    short_lipid_ids <- paste0("L", seq_along(original_lipid_names))
    lipid_name_map_diff <- data.frame(ID = short_lipid_ids, Name = original_lipid_names)

    # Save mapping so you can render it later
    output$lipid_name_lookup_diff <- DT::renderDataTable({
      lipid_name_map_diff
    }, options = list(pageLength = 10, scrollX = TRUE))

    # Rename columns in the heatmap matrix
    colnames(mat_comp) <- short_lipid_ids

    # Custom color palette with dark gray for NA values
    na_color <- "#2d2d2d"  # Dark gray/black for missing values

    # Get color palette based on user selection
    palette_choice <- input$heatmap_diff_color_palette
    if (is.null(palette_choice)) palette_choice <- "RdYlBu"

    if (palette_choice == "viridis") {
      main_colors <- viridisLite::viridis(255)
    } else if (palette_choice == "plasma") {
      main_colors <- viridisLite::plasma(255)
    } else {
      main_colors <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(255)
    }

    heatmaply(
      mat_comp,
      scale = "none",
      dendrogram = "none",
      showticklabels = c(TRUE, TRUE),
      colors = main_colors,
      na.value = na_color,
      xlab = "Lipids",
      ylab = "Condition",
      labRow = rownames(mat_comp),
      labCol = NULL,
      margins = c(60, 120, 40, 100)
    )
  })
  
  output$heatmap_summed <- renderPlotly({
    req(lipid_class_info, input$heatmap_summed_dataset, input$data_type_heatmap_summed)
    df <- get_dataset(input$heatmap_summed_dataset, input$data_type_heatmap_summed)

    # Extract only lipid cols with class info
    lipid_cols <- intersect(lipid_class_info$Lipids, colnames(df))
    class_map <- lipid_class_info %>% filter(Lipids %in% lipid_cols)
    
    # Keep only classes with multiple lipids
    class_counts <- table(class_map$Class)
    valid_classes <- names(class_counts[class_counts > 1])
    class_map <- class_map %>% filter(Class %in% valid_classes)
    
    df_sub <- df[, c("Compound_Name", class_map$Lipids)]
    
    # Get matrix of summed classes
    mat <- sapply(valid_classes, function(cls) {
      lipids <- class_map$Lipids[class_map$Class == cls]
      rowSums(df_sub[, lipids, drop = FALSE], na.rm = TRUE)
    })
    rownames(mat) <- df_sub$Compound_Name
    
    # Ranking by selected axis
    ranking_axis <- ifelse(input$heatmap_summed_mode == "Samples", 1, 2)
    diff_metric <- apply(mat, ranking_axis, function(x) {
      median_x <- median(x, na.rm = TRUE)
      mean(abs(x - median_x), na.rm = TRUE)
    })
    
    top_n <- input$top_n_summed
    selected <- names(sort(diff_metric, decreasing = TRUE))[1:min(top_n, length(diff_metric))]
    
    if (ranking_axis == 1) {
      mat <- mat[selected, , drop = FALSE]
    } else {
      mat <- mat[, selected, drop = FALSE]
    }
    
    mat_scaled <- scale(mat)
    
    heatmaply(
      mat_scaled,
      scale = "none",
      dendrogram = "both",
      colors = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(255),
      xlab = "Lipid Class",
      ylab = "Samples",
      labCol = colnames(mat_scaled),
      labRow = rownames(mat_scaled),
      margins = c(40, 120, 40, 120)
    )
  })
  
  
  
  
  
  ##======================================================================
  ##  TAB 7: GWAS
  ##======================================================================
  
  output$gwas_class_filter <- renderUI({
    req(input$gwas_dataset, input$gwas_threshold, input$gwas_trait_source)
    trait_meta <- get_gwas_trait_meta_table(
      input$gwas_dataset,
      source_choice = input$gwas_trait_source,
      p_threshold = input$gwas_threshold
    )
    classes <- sort(unique(unlist(trait_meta$ClassMembers, use.names = FALSE)))
    classes <- classes[classes != ""]
    choices <- c("All", classes)
    selectInput("gwas_class_choice", "Filter by Class:", choices = choices, selected = "All")
  })

  output$gwas_subclass_filter <- renderUI({
    req(input$gwas_dataset, input$gwas_threshold, input$gwas_trait_source, input$gwas_class_choice)
    trait_meta <- get_gwas_trait_meta_table(
      input$gwas_dataset,
      source_choice = input$gwas_trait_source,
      p_threshold = input$gwas_threshold
    )
    if (!is.null(input$gwas_class_choice) && input$gwas_class_choice != "All") {
      trait_meta <- trait_meta[vapply(trait_meta$ClassMembers, function(x) input$gwas_class_choice %in% x, logical(1)), , drop = FALSE]
    }
    subclasses <- sort(unique(unlist(trait_meta$SubclassMembers, use.names = FALSE)))
    choices <- subclasses[subclasses != ""]
    selectInput("gwas_subclass_choice", "Filter by Subclass:", choices = c("All", choices), selected = "All")
  })

  # Lipid dropdown based on dataset
  output$gwas_lipid_selector <- renderUI({
    req(input$gwas_dataset, input$gwas_threshold, input$gwas_trait_source, input$gwas_class_choice, input$gwas_subclass_choice)
    data_list <- get_annotation_data_list(input$gwas_dataset, p_threshold = input$gwas_threshold)
    trait_list <- names(data_list)
    if (length(trait_list) == 0) {
      return(tags$p("Add the GWAS annotation RDS files for this dataset to enable lipid selection."))
    }

    filtered_traits <- get_filtered_gwas_traits(
      dataset_name = input$gwas_dataset,
      source_choice = input$gwas_trait_source,
      class_choice = input$gwas_class_choice,
      subclass_choice = input$gwas_subclass_choice,
      p_threshold = input$gwas_threshold
    )
    
    hit_counts <- vapply(data_list, function(df) if (is.null(df)) 0L else nrow(df), integer(1))
    ordered_lipids <- c(
      filtered_traits[hit_counts[filtered_traits] > 0],
      filtered_traits[hit_counts[filtered_traits] == 0]
    )
    ordered_lipids <- unique(ordered_lipids)
    if (length(ordered_lipids) == 0) {
      return(tags$p("No GWAS traits match the selected class/subclass filter."))
    }

    lipid_choices <- c(
      setNames("__ALL__", paste0("All (", length(ordered_lipids), " traits)")),
      setNames(ordered_lipids, paste0(ordered_lipids, " (", hit_counts[ordered_lipids], " hits)"))
    )
    
    selectInput("selected_lipid", "Select Lipid:", choices = lipid_choices, selected = "__ALL__")
  })

  selected_gwas_trait_meta <- reactive({
    req(input$gwas_dataset, input$gwas_threshold, input$gwas_trait_source, input$selected_lipid)
    if (input$selected_lipid == "__ALL__") {
      return(NULL)
    }
    hit <- get_gwas_trait_meta_table(
      input$gwas_dataset,
      trait_names = input$selected_lipid,
      source_choice = input$gwas_trait_source,
      p_threshold = input$gwas_threshold
    )
    if (nrow(hit) == 0) {
      return(NULL)
    }
    hit[1, c("Class", "Subclass"), drop = FALSE]
  })

  output$gwas_trait_class_info <- renderUI({
    req(input$gwas_dataset, input$gwas_threshold, input$gwas_trait_source, input$selected_lipid)
    if (input$selected_lipid == "__ALL__") {
      traits <- get_filtered_gwas_traits(
        dataset_name = input$gwas_dataset,
        source_choice = input$gwas_trait_source,
        class_choice = input$gwas_class_choice,
        subclass_choice = input$gwas_subclass_choice,
        p_threshold = input$gwas_threshold
      )
      source_label <- dplyr::case_when(
        input$gwas_trait_source == "all" ~ "All",
        input$gwas_trait_source == "individual" ~ "Individual",
        input$gwas_trait_source == "sum_ratio" ~ "Sum/Ratio",
        TRUE ~ input$gwas_trait_source
      )
      return(tags$p(
        style = "margin-bottom: 10px; color: #555;",
        tags$strong("Filtered traits: "), length(traits), "  |  ",
        tags$strong("Source: "), source_label, "  |  ",
        tags$strong("Threshold: "), paste0("-log10(p) >= ", input$gwas_threshold)
      ))
    }

    meta <- selected_gwas_trait_meta()
    if (is.null(meta)) {
      return(NULL)
    }

    class_text <- ifelse(meta$Class[[1]] == "", "NA", meta$Class[[1]])
    subclass_text <- ifelse(meta$Subclass[[1]] == "", "NA", meta$Subclass[[1]])

    tags$p(
      style = "margin-bottom: 10px; color: #555;",
      tags$strong("Class: "), class_text, "  |  ",
      tags$strong("Subclass: "), subclass_text
    )
  })
  
  # Manhattan plot (dynamic path based on dataset)
  output$manhattan_image <- renderUI({
    req(input$gwas_dataset, input$selected_lipid)
    if (input$selected_lipid == "__ALL__") {
      return(tags$p("Select a specific lipid to view Manhattan plot."))
    }
    lipid <- input$selected_lipid
    
    # Build the appropriate directory path
    dir_name <- paste0("manhattans/", input$gwas_dataset)
    img_path <- file.path(dir_name, paste0(lipid, ".jpg"))
    
    # Check if file exists inside www/
    if (!file.exists(file.path("www", img_path))) {
      return(tags$p("No Manhattan plot available for this lipid."))
    }
    
    tags$img(
      src = img_path,
      width = "100%",
      style = "border:1px solid #ddd; padding:10px;"
    )
  })
  
  
  # GWAS Table
  output$gwas_table <- renderDT({
    req(input$gwas_dataset, input$gwas_threshold, input$gwas_trait_source, input$selected_lipid)
    data_list <- get_annotation_data_list(input$gwas_dataset, p_threshold = input$gwas_threshold)

    if (input$selected_lipid == "__ALL__") {
      traits <- get_filtered_gwas_traits(
        dataset_name = input$gwas_dataset,
        source_choice = input$gwas_trait_source,
        class_choice = input$gwas_class_choice,
        subclass_choice = input$gwas_subclass_choice,
        p_threshold = input$gwas_threshold
      )
      rows <- lapply(traits, function(trait) {
        df <- data_list[[trait]]
        if (is.null(df) || nrow(df) == 0) {
          return(NULL)
        }
        df <- df %>% mutate(Trait = trait, .before = GeneID)
        meta <- get_gwas_trait_meta_table(
          input$gwas_dataset,
          trait_names = trait,
          source_choice = input$gwas_trait_source,
          p_threshold = input$gwas_threshold
        )
        df %>% mutate(
          Class = meta$Class[[1]],
          Subclass = meta$Subclass[[1]],
          .before = GeneID
        )
      })
      rows <- rows[!vapply(rows, is.null, logical(1))]
      validate(need(length(rows) > 0, "No genes available for the selected class/subclass filter."))
      df_annotated <- bind_rows(rows) %>%
        left_join(gene_annotation, by = "GeneID")
      return(datatable(df_annotated, options = list(scrollX = TRUE, pageLength = 10)))
    }

    df <- data_list[[input$selected_lipid]]
    validate(need(!is.null(df), "No annotation table found for this lipid."))
    validate(need(nrow(df) > 0, "No genes passed the GWAS cutoff for this lipid in the loaded RDS."))

    df_annotated <- df %>%
      left_join(gene_annotation, by = "GeneID")

    meta <- selected_gwas_trait_meta()
    if (!is.null(meta)) {
      df_annotated <- df_annotated %>%
        mutate(
          Class = meta$Class[[1]],
          Subclass = meta$Subclass[[1]],
          .before = GeneID
        )
    }

    datatable(df_annotated, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  ##======================================================================
  ##  TAB 8: Gene Hits
  ##======================================================================
  
  # Cascading filters for trait source -> class -> subclass -> lipid
  output$hit_class_filter <- renderUI({
    req(input$hit_dataset, input$hit_threshold, input$hit_trait_source)
    trait_meta <- get_gwas_trait_meta_table(
      input$hit_dataset,
      source_choice = input$hit_trait_source,
      p_threshold = input$hit_threshold
    )
    classes <- sort(unique(unlist(trait_meta$ClassMembers, use.names = FALSE)))
    classes <- classes[classes != ""]
    selectInput("hit_class_choice", "Filter by Class:", choices = c("All", classes), selected = "All")
  })

  output$hit_subclass_filter <- renderUI({
    req(input$hit_dataset, input$hit_threshold, input$hit_trait_source, input$hit_class_choice)
    trait_meta <- get_gwas_trait_meta_table(
      input$hit_dataset,
      source_choice = input$hit_trait_source,
      p_threshold = input$hit_threshold
    )
    if (!is.null(input$hit_class_choice) && input$hit_class_choice != "All") {
      trait_meta <- trait_meta[vapply(trait_meta$ClassMembers, function(x) input$hit_class_choice %in% x, logical(1)), , drop = FALSE]
    }
    subclasses <- sort(unique(unlist(trait_meta$SubclassMembers, use.names = FALSE)))
    subclasses <- subclasses[subclasses != ""]
    selectInput("hit_subclass_choice", "Filter by Subclass:", choices = c("All", subclasses), selected = "All")
  })

  output$hit_lipid_filter <- renderUI({
    req(input$hit_dataset, input$hit_threshold, input$hit_trait_source, input$hit_class_choice, input$hit_subclass_choice)
    data_list <- get_annotation_data_list(input$hit_dataset, p_threshold = input$hit_threshold)
    validate(need(length(data_list) > 0, "No GWAS annotations loaded for the selected dataset."))

    traits <- get_filtered_gwas_traits(
      dataset_name = input$hit_dataset,
      source_choice = input$hit_trait_source,
      class_choice = input$hit_class_choice,
      subclass_choice = input$hit_subclass_choice,
      p_threshold = input$hit_threshold
    )

    if (length(traits) == 0) {
      return(selectInput("hit_lipid_choice", "Filter by Lipid:", choices = "All", selected = "All"))
    }

    hit_counts <- vapply(data_list[traits], function(df) if (is.null(df)) 0L else nrow(df), integer(1))
    trait_choices <- c(
      "All" = "All",
      setNames(traits, paste0(traits, " (", hit_counts[traits], " hits)"))
    )
    selectInput("hit_lipid_choice", "Filter by Lipid:", choices = trait_choices, selected = "All")
  })

  output$gene_hit_table <- renderDT({
    req(
      input$hit_dataset,
      input$hit_threshold,
      input$hit_trait_source,
      input$hit_class_choice,
      input$hit_subclass_choice,
      input$hit_lipid_choice
    )
    data_list <- get_annotation_data_list(input$hit_dataset, p_threshold = input$hit_threshold)
    validate(need(length(data_list) > 0, "No GWAS annotations loaded for the selected dataset."))
    threshold <- as.numeric(input$hit_threshold)

    trait_subset <- get_filtered_gwas_traits(
      dataset_name = input$hit_dataset,
      source_choice = input$hit_trait_source,
      class_choice = input$hit_class_choice,
      subclass_choice = input$hit_subclass_choice,
      lipid_choice = input$hit_lipid_choice,
      p_threshold = input$hit_threshold
    )

    # Combine hits for selected trait subset
    hit_rows <- lapply(trait_subset, function(trait) {
      df <- data_list[[trait]]
      if (!is.null(df) && all(c("GeneID", "log(p)") %in% colnames(df))) {
        out <- df[df$`log(p)` >= threshold, c("GeneID", "log(p)"), drop = FALSE]
        out$Trait <- trait
        out
      }
    })
    hit_rows <- hit_rows[!vapply(hit_rows, is.null, logical(1))]
    validate(need(length(hit_rows) > 0, "No genes passed the selected threshold/filter combination."))

    combined_df <- bind_rows(hit_rows)

    gene_summary <- combined_df %>%
      group_by(GeneID) %>%
      summarise(
        Occurrences = n(),
        Highest_logp = round(max(`log(p)`, na.rm = TRUE), 2),
        .groups = "drop"
      ) %>%
      arrange(desc(Occurrences), desc(Highest_logp)) %>%
      left_join(gene_annotation, by = "GeneID")
    
    datatable(gene_summary, options = list(scrollX = TRUE, pageLength = 10))
  })

  output$phenotype_hit_table <- renderDT({
    req(
      input$hit_dataset,
      input$hit_threshold,
      input$hit_trait_source,
      input$hit_class_choice,
      input$hit_subclass_choice,
      input$hit_lipid_choice
    )
    data_list <- get_annotation_data_list(input$hit_dataset, p_threshold = input$hit_threshold)
    validate(need(length(data_list) > 0, "No GWAS annotations loaded for the selected dataset."))
    threshold <- as.numeric(input$hit_threshold)

    trait_subset <- get_filtered_gwas_traits(
      dataset_name = input$hit_dataset,
      source_choice = input$hit_trait_source,
      class_choice = input$hit_class_choice,
      subclass_choice = input$hit_subclass_choice,
      lipid_choice = input$hit_lipid_choice,
      p_threshold = input$hit_threshold
    )
    trait_meta <- get_gwas_trait_meta_table(
      input$hit_dataset,
      trait_names = trait_subset,
      source_choice = input$hit_trait_source,
      p_threshold = input$hit_threshold
    ) %>%
      transmute(Trait = Trait, Class = Class, Subclass = Subclass)
    source_map <- get_annotation_source_map(input$hit_dataset, p_threshold = input$hit_threshold)

    summary_rows <- lapply(trait_subset, function(trait) {
      df <- data_list[[trait]]
      if (is.null(df) || !all(c("GeneID", "log(p)") %in% colnames(df))) {
        return(NULL)
      }
      df_sig <- df[df$`log(p)` >= threshold, c("GeneID", "log(p)"), drop = FALSE]
      if (nrow(df_sig) == 0) {
        return(NULL)
      }
      idx <- which.max(df_sig$`log(p)`)
      trait_info <- trait_meta %>% filter(Trait == trait)
      data.frame(
        Phenotype = trait,
        TraitSource = if (!is.null(source_map)) source_map[[trait]] else "unknown",
        Class = if (nrow(trait_info) > 0) trait_info$Class[[1]] else "",
        Subclass = if (nrow(trait_info) > 0) trait_info$Subclass[[1]] else "",
        SignificantGenes = nrow(df_sig),
        Highest_logp = round(max(df_sig$`log(p)`, na.rm = TRUE), 2),
        TopGeneID = df_sig$GeneID[[idx]],
        stringsAsFactors = FALSE
      )
    })
    summary_rows <- summary_rows[!vapply(summary_rows, is.null, logical(1))]
    validate(need(length(summary_rows) > 0, "No phenotypes passed the selected threshold/filter combination."))

    summary_df <- bind_rows(summary_rows) %>%
      arrange(desc(Highest_logp), desc(SignificantGenes), Phenotype)

    datatable(summary_df, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  
  ##======================================================================
  ##  TAB 9: Correlations
  ##======================================================================
  
  # output$lipid1_selector <- renderUI({
  #   df <- dataset_list[[input$cor_dataset]]
  #   lipid_names <- colnames(df)[-1]
  #   selectInput("lipid1", "Select Lipid 1:", choices = lipid_names)
  # })
  # 
  # output$lipid2_selector <- renderUI({
  #   df <- dataset_list[[input$cor_dataset]]
  #   lipid_names <- colnames(df)[-1]
  #   selectInput("lipid2", "Select Lipid 2:", choices = lipid_names)
  # })
  # 
  # output$correlation_plot <- renderPlot({
  #   req(input$cor_dataset, input$lipid1, input$lipid2)
  #   
  #   df <- dataset_list[[input$cor_dataset]]
  #   
  #   validate(
  #     need(all(c(input$lipid1, input$lipid2) %in% colnames(df)), "Selected lipids not found in dataset")
  #   )
  #   
  #   # Clean: filter out rows where either value is NA
  #   df_clean <- df %>%
  #     dplyr::filter(!is.na(.data[[input$lipid1]]) & !is.na(.data[[input$lipid2]]))
  #   
  #   # Dynamically assign variables for plotting
  #   ggplot(df_clean, aes(x = .data[[input$lipid1]], y = .data[[input$lipid2]])) +
  #     geom_point(alpha = 0.7) +
  #     geom_smooth(method = "lm", color = "blue", se = FALSE) +
  #     labs(
  #       title = paste("Correlation between", input$lipid1, "and", input$lipid2),
  #       x = input$lipid1,
  #       y = input$lipid2
  #     ) +
  #     theme_minimal() +
  #     theme_minimal(base_size = 20) +
  #     theme(
  #       panel.grid.major = element_blank(),  # Remove major grid lines
  #       panel.grid.minor = element_blank(),  # Remove minor grid lines
  #       axis.line = element_line(),          # Add base axes
  #       axis.text.x = element_text(angle = 45, hjust = 1),
  #       legend.position = "none"
  #     )
  # })
  # 
  # # UI for lipid choices in the selected dataset
  # output$geno_lipid_ui <- renderUI({
  #   req(input$geno_dataset)
  #   df <- dataset_list[[input$geno_dataset]]
  #   lipid_choices <- colnames(df)[-1]
  #   selectInput("geno_lipid", "Choose Lipid:", choices = lipid_choices)
  # })
  # 
  # # UI for slider based on selected lipid range
  # output$geno_slider_ui <- renderUI({
  #   req(input$geno_dataset, input$geno_lipid)
  #   df <- dataset_list[[input$geno_dataset]]
  #   vals <- df[[input$geno_lipid]]
  #   sliderInput("geno_range", "Select Lipid Value Range:",
  #               min = floor(min(vals, na.rm = TRUE)),
  #               max = ceiling(max(vals, na.rm = TRUE)),
  #               value = c(floor(min(vals, na.rm = TRUE)), ceiling(max(vals, na.rm = TRUE))),
  #               step = 0.01)
  # })
  
  output$lipid1_selector <- renderUI({
    req(input$cor_dataset, input$data_type_cor)
    df <- get_dataset(input$cor_dataset, input$data_type_cor)
    req(df)
    lipids <- colnames(df)[-1]

    # Filter by class if selected
    class_filter <- if (!is.null(input$cor_class_filter)) input$cor_class_filter else "all"
    if (class_filter != "all") {
      lipid_classes <- sapply(lipids, get_lipid_class)
      lipids <- lipids[lipid_classes == class_filter]
    }

    if (length(lipids) == 0) lipids <- colnames(df)[-1]  # fallback
    selectInput("lipid1", "Select Lipid 1:", choices = lipids, selected = lipids[1])
  })

  output$lipid2_selector <- renderUI({
    req(input$cor_dataset, input$data_type_cor)
    df <- get_dataset(input$cor_dataset, input$data_type_cor)
    req(df)
    lipids <- colnames(df)[-1]

    # Filter by class if selected
    class_filter <- if (!is.null(input$cor_class_filter)) input$cor_class_filter else "all"
    if (class_filter != "all") {
      lipid_classes <- sapply(lipids, get_lipid_class)
      lipids <- lipids[lipid_classes == class_filter]
    }

    if (length(lipids) == 0) lipids <- colnames(df)[-1]  # fallback
    selectInput("lipid2", "Select Lipid 2:", choices = lipids, selected = lipids[min(2, length(lipids))])
  })

  output$correlation_plot <- renderPlot({
    req(input$cor_dataset, input$lipid1, input$lipid2, input$data_type_cor)
    df <- get_dataset(input$cor_dataset, input$data_type_cor)
    validate(need(all(c(input$lipid1, input$lipid2) %in% colnames(df)), "Selected lipids not found"))

    df_clean <- df %>% dplyr::filter(!is.na(.data[[input$lipid1]]) & !is.na(.data[[input$lipid2]]))

    # Dynamic axis label based on data type
    data_label <- paste0(" (", data_type_label(input$data_type_cor), ")")

    p <- ggplot(df_clean, aes(x = .data[[input$lipid1]], y = .data[[input$lipid2]])) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", color = "blue", se = FALSE) +
      labs(title = paste("Correlation between", input$lipid1, "and", input$lipid2, data_label),
           x = input$lipid1, y = input$lipid2) +
      plot_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    plot_store$correlation_plot <- p
    p
  })
  
  output$class_correlation_heatmap <- renderPlotly({
    req(input$class_cor_dataset, input$data_type_class_cor)
    df <- get_dataset(input$class_cor_dataset, input$data_type_class_cor)
    
    lipid_classes <- lipid_class_info %>%
      dplyr::filter(Lipids %in% colnames(df)) %>%
      group_by(Class) %>%
      summarise(Lipids = list(Lipids), .groups = 'drop')
    
    class_matrix <- sapply(lipid_classes$Lipids, function(lip_set) {
      rowSums(df[, intersect(unlist(lip_set), colnames(df)), drop = FALSE], na.rm = TRUE)
    })
    
    class_matrix <- as.data.frame(class_matrix)
    colnames(class_matrix) <- lipid_classes$Class
    rownames(class_matrix) <- df[[1]]
    
    validate(need(ncol(class_matrix) > 1, "Not enough lipid classes for correlation."))
    
    corr_matrix <- cor(class_matrix, use = "pairwise.complete.obs", method = "pearson")
    
    heatmaply::heatmaply(
      corr_matrix,
      xlab = "Lipid Classes",
      ylab = "Lipid Classes",
      main = "Summed Class Correlation",
      colors = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100)
    )
  })
  
  output$individual_correlation_heatmap <- renderPlotly({
    req(input$corheat_dataset, input$top_mad_lipids, input$data_type_corheat)

    df <- get_dataset(input$corheat_dataset, input$data_type_corheat)
    mat <- df[, -1]
    rownames(mat) <- df[[1]]
    
    # Compute MAD per lipid
    mad_vals <- apply(mat, 2, function(x) mean(abs(x - median(x, na.rm = TRUE)), na.rm = TRUE))
    
    top_lipids <- names(sort(mad_vals, decreasing = TRUE))[1:min(input$top_mad_lipids, length(mad_vals))]
    mat_top <- mat[, top_lipids, drop = FALSE]
    
    # Compute correlation matrix
    corr_mat <- cor(mat_top, use = "pairwise.complete.obs", method = "pearson")
    
    validate(need(ncol(corr_mat) > 1, "Not enough lipids for correlation heatmap."))
    
    heatmaply::heatmaply(
      corr_mat,
      xlab = "Lipids", ylab = "Lipids", main = "Top Correlated Lipids (Pearson)",
      colors = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(255),
      margins = c(60, 100, 40, 100)
    )
  })
  
  
  
  ##======================================================================
  ##  TAB 10: Select Genotype
  ##======================================================================
  
  # UI: Dataset and Lipid Selectors
  output$geno_dataset_ui <- renderUI({
    selectInput("geno_dataset", "Choose Dataset:", choices = names(dataset_list))
  })

  # Second dataset selector for Genotype Viewer tab
  output$geno_dataset_ui2 <- renderUI({
    selectInput("geno_dataset2", "Choose Dataset:", choices = names(dataset_list))
  })

  output$geno_lipid_ui <- renderUI({
    req(input$geno_dataset, input$data_type_geno)
    df <- get_dataset(input$geno_dataset, input$data_type_geno)
    all_lipids <- colnames(df)[-1]
    class_options <- paste0("SUM_", sort(unique(lipid_class_info$Class)))
    selectInput("geno_lipid", "Choose Lipid:", choices = c(class_options, all_lipids))
  })

  # UI: Slider for Custom Range
  output$geno_slider_ui <- renderUI({
    req(input$geno_dataset, input$geno_lipid, input$data_type_geno)
    df <- get_dataset(input$geno_dataset, input$data_type_geno)
    
    values <- if (startsWith(input$geno_lipid, "SUM_")) {
      class_name <- sub("SUM_", "", input$geno_lipid)
      lipids <- lipid_class_info$Lipids[lipid_class_info$Class == class_name]
      valid <- intersect(lipids, colnames(df))
      rowSums(df[, valid, drop = FALSE], na.rm = TRUE)
    } else {
      df[[input$geno_lipid]]
    }
    
    sliderInput("geno_range", "Select Lipid Value Range:",
                min = floor(min(values, na.rm = TRUE)),
                max = ceiling(max(values, na.rm = TRUE)),
                value = quantile(values, c(0.25, 0.75), na.rm = TRUE))
  })
  
  # Logic: Handle Filter Button
  observeEvent(input$filter_genotypes, {
    req(input$geno_dataset, input$geno_lipid, input$data_type_geno)
    df <- get_dataset(input$geno_dataset, input$data_type_geno)
    
    values <- if (startsWith(input$geno_lipid, "SUM_")) {
      class_name <- sub("SUM_", "", input$geno_lipid)
      lipids <- lipid_class_info$Lipids[lipid_class_info$Class == class_name]
      valid <- intersect(lipids, colnames(df))
      rowSums(df[, valid, drop = FALSE], na.rm = TRUE)
    } else {
      df[[input$geno_lipid]]
    }
    
    selected_genos <- NULL
    
    if (input$geno_filter_method == "Lowest Quartile") {
      q1 <- quantile(values, 0.25, na.rm = TRUE)
      selected_genos <- df[[1]][values <= q1]
      
    } else if (input$geno_filter_method == "Highest Quartile") {
      q3 <- quantile(values, 0.75, na.rm = TRUE)
      selected_genos <- df[[1]][values >= q3]
      
    } else {
      req(input$geno_range)
      selected_genos <- df[[1]][values >= input$geno_range[1] & values <= input$geno_range[2]]
    }
    
    output$genotype_output <- renderPrint({
      if (length(selected_genos) == 0) {
        "No genotypes found for the selected filter."
      } else {
        cat("Genotypes in selected filter:\n\n")
        print(selected_genos)
      }
    })
  })
  
  ##======================================================================
  ##  TAB 10: Genotype Viewer
  ##======================================================================
  
  # UI: Genotype selector and lipid selector
  output$geno_viewer_ui <- renderUI({
    req(input$geno_dataset2, input$data_type_geno_view)
    df <- get_dataset(input$geno_dataset2, input$data_type_geno_view)
    selectInput("geno_view_id", "Select Genotype:", choices = df[[1]])
  })

  output$geno_view_lipid_ui <- renderUI({
    req(input$geno_dataset2, input$data_type_geno_view)
    df <- get_dataset(input$geno_dataset2, input$data_type_geno_view)
    all_lipids <- colnames(df)[-1]
    class_options <- paste0("SUM_", sort(unique(lipid_class_info$Class)))
    selectInput("geno_view_lipid", "Select Lipid:", choices = c(class_options, all_lipids))
  })

  # Plot viewer (barplot for both conditions)
  output$geno_view_plot <- renderPlot({
    req(input$geno_view_id, input$geno_view_lipid, input$data_type_geno_view)

    get_lipid_value <- function(df, genotype, lipid) {
      if (startsWith(lipid, "SUM_")) {
        class_name <- sub("SUM_", "", lipid)
        lipids <- lipid_class_info$Lipids[lipid_class_info$Class == class_name]
        valid <- intersect(lipids, colnames(df))
        rowSums(df[df[[1]] == genotype, valid, drop = FALSE], na.rm = TRUE)
      } else {
        df[df[[1]] == genotype, lipid, drop = TRUE]
      }
    }

    # Get all values for the selected lipid across all genotypes (for reference lines)
    get_all_lipid_values <- function(df, lipid) {
      if (startsWith(lipid, "SUM_")) {
        class_name <- sub("SUM_", "", lipid)
        lipids <- lipid_class_info$Lipids[lipid_class_info$Class == class_name]
        valid <- intersect(lipids, colnames(df))
        rowSums(df[, valid, drop = FALSE], na.rm = TRUE)
      } else {
        as.numeric(df[[lipid]])
      }
    }

    values <- sapply(names(dataset_list), function(dataset) {
      get_lipid_value(get_dataset(dataset, input$data_type_geno_view), input$geno_view_id, input$geno_view_lipid)
    })

    df_bar <- data.frame(Condition = names(values), Value = as.numeric(values))

    # Dynamic y-axis label based on data type
    y_label <- if (input$data_type_geno_view == "raw") "Raw Signal Intensity" else data_type_label(input$data_type_geno_view)
    data_label <- paste0(" (", data_type_label(input$data_type_geno_view), ")")

    p <- ggplot(df_bar, aes(x = Condition, y = Value, fill = Condition)) +
      geom_col(width = 0.6) +
      scale_fill_manual(values = c("control" = "#FF7F0E", "lowinput" = "#2CA02C")) +
      labs(title = paste("Lipid Level for", input$geno_view_id, data_label), y = y_label, x = "Condition") +
      plot_theme +
      theme(legend.position = "none")

    # Add reference lines based on selection
    ref_lines <- if (!is.null(input$geno_ref_lines)) input$geno_ref_lines else "none"

    if (ref_lines != "none") {
      # Calculate statistics for each condition
      ref_data <- lapply(names(dataset_list), function(dataset) {
        all_vals <- get_all_lipid_values(get_dataset(dataset, input$data_type_geno_view), input$geno_view_lipid)
        all_vals <- all_vals[!is.na(all_vals)]
        list(
          condition = dataset,
          mean = mean(all_vals, na.rm = TRUE),
          median = median(all_vals, na.rm = TRUE),
          q1 = quantile(all_vals, 0.25, na.rm = TRUE),
          q3 = quantile(all_vals, 0.75, na.rm = TRUE)
        )
      })

      ref_df <- do.call(rbind, lapply(ref_data, function(x) {
        data.frame(
          Condition = x$condition,
          mean = x$mean,
          median = x$median,
          q1 = x$q1,
          q3 = x$q3
        )
      }))

      # Add reference lines based on selection
      if (ref_lines == "mean") {
        p <- p + geom_errorbar(data = ref_df, aes(x = Condition, ymin = mean, ymax = mean),
                               width = 0.4, linewidth = 1, color = "red", linetype = "dashed", inherit.aes = FALSE) +
          geom_text(data = ref_df, aes(x = Condition, y = mean, label = "Mean"),
                    vjust = -0.5, color = "red", size = 3, inherit.aes = FALSE)
      } else if (ref_lines == "median") {
        p <- p + geom_errorbar(data = ref_df, aes(x = Condition, ymin = median, ymax = median),
                               width = 0.4, linewidth = 1, color = "blue", linetype = "dashed", inherit.aes = FALSE) +
          geom_text(data = ref_df, aes(x = Condition, y = median, label = "Median"),
                    vjust = -0.5, color = "blue", size = 3, inherit.aes = FALSE)
      } else if (ref_lines == "q1") {
        p <- p + geom_errorbar(data = ref_df, aes(x = Condition, ymin = q1, ymax = q1),
                               width = 0.4, linewidth = 1, color = "purple", linetype = "dashed", inherit.aes = FALSE) +
          geom_text(data = ref_df, aes(x = Condition, y = q1, label = "Q1"),
                    vjust = -0.5, color = "purple", size = 3, inherit.aes = FALSE)
      } else if (ref_lines == "q3") {
        p <- p + geom_errorbar(data = ref_df, aes(x = Condition, ymin = q3, ymax = q3),
                               width = 0.4, linewidth = 1, color = "darkgreen", linetype = "dashed", inherit.aes = FALSE) +
          geom_text(data = ref_df, aes(x = Condition, y = q3, label = "Q3"),
                    vjust = -0.5, color = "darkgreen", size = 3, inherit.aes = FALSE)
      } else if (ref_lines == "all") {
        p <- p +
          geom_errorbar(data = ref_df, aes(x = Condition, ymin = q1, ymax = q1),
                        width = 0.4, linewidth = 0.8, color = "purple", linetype = "dotted", inherit.aes = FALSE) +
          geom_errorbar(data = ref_df, aes(x = Condition, ymin = median, ymax = median),
                        width = 0.4, linewidth = 1, color = "blue", linetype = "dashed", inherit.aes = FALSE) +
          geom_errorbar(data = ref_df, aes(x = Condition, ymin = mean, ymax = mean),
                        width = 0.4, linewidth = 1, color = "red", linetype = "dashed", inherit.aes = FALSE) +
          geom_errorbar(data = ref_df, aes(x = Condition, ymin = q3, ymax = q3),
                        width = 0.4, linewidth = 0.8, color = "darkgreen", linetype = "dotted", inherit.aes = FALSE) +
          geom_text(data = ref_df, aes(x = Condition, y = q1, label = "Q1"),
                    vjust = 1.5, color = "purple", size = 2.5, inherit.aes = FALSE) +
          geom_text(data = ref_df, aes(x = Condition, y = median, label = "Med"),
                    vjust = -0.5, color = "blue", size = 2.5, inherit.aes = FALSE) +
          geom_text(data = ref_df, aes(x = Condition, y = mean, label = "Mean"),
                    vjust = 1.5, color = "red", size = 2.5, inherit.aes = FALSE) +
          geom_text(data = ref_df, aes(x = Condition, y = q3, label = "Q3"),
                    vjust = -0.5, color = "darkgreen", size = 2.5, inherit.aes = FALSE)
      }
    }

    plot_store$geno_view_plot <- p
    p
  })
  
  ##======================================================================
  ##  TAB 11: Venn Diagram
  ##======================================================================
  
  # UI ----
  output$venn_lipid_type_ui <- renderUI({
    selectInput("venn_type", "Lipid Type:", choices = c("Individual", "Class"))
  })
  
  output$venn_class_selector <- renderUI({
    req(input$venn_type == "Class")
    selectInput("venn_class_choice", "Select Lipid Class:", 
                choices = sort(unique(lipid_class_info$Class)))
  })
  
  # Server Logic ----
  output$venn_plot <- renderPlot({
    req(input$venn_type)

    df_control <- dataset_list$control
    df_lowinput <- dataset_list$lowinput
    
    # Get lipid sets based on type
    if (input$venn_type == "Individual") {
      set_control <- colnames(df_control)[-1]  # skip Sample column
      set_lowinput <- colnames(df_lowinput)[-1]
    } else {
      req(input$venn_class_choice)
      class_lipids <- lipid_class_info$Lipids[lipid_class_info$Class == input$venn_class_choice]
      
      set_control <- intersect(class_lipids, colnames(df_control))
      set_lowinput <- intersect(class_lipids, colnames(df_lowinput))
    }
    
    venn.plot <- draw.pairwise.venn(
      area1 = length(set_control),
      area2 = length(set_lowinput),
      cross.area = length(intersect(set_control, set_lowinput)),
      category = c("Control", "Low Input"),
      fill = c("skyblue", "orange"),
      alpha = c(0.5, 0.5),
      lty = "blank",
      cex = 2,
      cat.cex = 1.5,
      cat.pos = c(-20, 20),
      cat.dist = 0.05,
      scaled = TRUE
    )
    
    grid.newpage()
    grid.draw(venn.plot)
  })
  
  output$venn_table <- renderDT({
    req(input$venn_type)
    
    df_control <- dataset_list$control
    df_lowinput <- dataset_list$lowinput
    
    if (input$venn_type == "Individual") {
      set_control <- colnames(df_control)[-1]
      set_lowinput <- colnames(df_lowinput)[-1]
    } else {
      req(input$venn_class_choice)
      class_lipids <- lipid_class_info$Lipids[lipid_class_info$Class == input$venn_class_choice]
      
      set_control <- intersect(class_lipids, colnames(df_control))
      set_lowinput <- intersect(class_lipids, colnames(df_lowinput))
    }
    
    only_control <- setdiff(set_control, set_lowinput)
    only_lowinput <- setdiff(set_lowinput, set_control)
    shared <- intersect(set_control, set_lowinput)
    
    max_len <- max(length(only_control), length(shared), length(only_lowinput))
    
    venn_df <- data.frame(
      `Only Control` = c(only_control, rep(NA, max_len - length(only_control))),
      Shared = c(shared, rep(NA, max_len - length(shared))),
      `Only Low Input` = c(only_lowinput, rep(NA, max_len - length(only_lowinput)))
    )
    
    datatable(venn_df, options = list(pageLength = 10, scrollX = TRUE))
  })

  ##======================================================================
  ##  PCA Variance Information
  ##======================================================================
  output$pca_variance_info <- renderUI({
    req(input$pca_dataset, input$data_type_pca)
    pca_data_mode <- if (input$data_type_pca == "clr") "tic" else input$data_type_pca
    df <- get_dataset(input$pca_dataset, pca_data_mode)
    pca_bundle <- compute_pca_bundle(
      df = df,
      mode = input$pca_mode,
      data_type = input$data_type_pca,
      class_filter = if (input$pca_mode == "Lipids") input$pca_class_filter else "All",
      center = input$pca_center,
      scale = input$pca_scale
    )
    if (is.null(pca_bundle) || is.null(pca_bundle$pve)) {
      return(tags$p(style = "color: #666; margin-top: 10px;", "Variance information not available for this selection."))
    }

    var_explained <- 100 * pca_bundle$pve[1:min(5, length(pca_bundle$pve))]

    tags$div(
      style = "background: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 10px;",
      tags$strong("Variance Explained:"),
      tags$ul(
        lapply(seq_along(var_explained), function(i) {
          tags$li(paste0("PC", i, ": ", round(var_explained[i], 1), "%"))
        })
      ),
      tags$p(style = "color: #666; font-size: 12px;",
             paste0("Total (first 5 PCs): ", round(sum(var_explained), 1), "%")),
      tags$p(style = "color: #666; font-size: 12px; margin-bottom: 0;",
             paste0("Method: ", pca_bundle$data_label))
    )
  })

  ##======================================================================
  ##  VOLCANO PLOT
  ##======================================================================

  # Reactive: compute differential analysis
  volcano_data <- reactive({
    req(input$volcano_data_type)

    df_control <- get_dataset("control", input$volcano_data_type)
    df_lowinput <- get_dataset("lowinput", input$volcano_data_type)

    # Get common lipids
    common_lipids <- intersect(colnames(df_control)[-1], colnames(df_lowinput)[-1])

    # Calculate fold change and p-value for each lipid
    results <- lapply(common_lipids, function(lipid) {
      ctrl_vals <- as.numeric(df_control[[lipid]])
      low_vals <- as.numeric(df_lowinput[[lipid]])

      # Remove NAs
      ctrl_vals <- ctrl_vals[!is.na(ctrl_vals)]
      low_vals <- low_vals[!is.na(low_vals)]

      if (length(ctrl_vals) < 3 || length(low_vals) < 3) {
        return(NULL)
      }

      # Mean values
      mean_ctrl <- mean(ctrl_vals, na.rm = TRUE)
      mean_low <- mean(low_vals, na.rm = TRUE)

      # Avoid log of zero
      if (mean_ctrl <= 0) mean_ctrl <- 1e-10
      if (mean_low <= 0) mean_low <- 1e-10

      # Log2 fold change (low input vs control)
      log2fc <- log2(mean_low / mean_ctrl)

      # Wilcoxon test (more robust than t-test)
      pval <- tryCatch({
        test_result <- wilcox.test(low_vals, ctrl_vals)
        test_result$p.value
      }, error = function(e) {
        1
      })

      data.frame(
        Lipid = lipid,
        Log2FC = log2fc,
        PValue = pval,
        Mean_Control = mean_ctrl,
        Mean_LowInput = mean_low,
        stringsAsFactors = FALSE
      )
    })

    results_df <- do.call(rbind, results[!sapply(results, is.null)])
    results_df$NegLog10P <- -log10(results_df$PValue)
    results_df$NegLog10P[is.infinite(results_df$NegLog10P)] <- max(results_df$NegLog10P[is.finite(results_df$NegLog10P)], na.rm = TRUE) + 1

    # Add lipid class
    results_df$Class <- sapply(results_df$Lipid, get_lipid_class)

    # Determine significance
    fc_thresh <- input$volcano_fc_threshold
    pval_thresh <- input$volcano_pval_threshold

    results_df$Significant <- ifelse(
      abs(results_df$Log2FC) >= fc_thresh & results_df$NegLog10P >= pval_thresh,
      ifelse(results_df$Log2FC > 0, "Up in Low Input", "Down in Low Input"),
      "Not Significant"
    )

    results_df
  })

  output$volcano_plot <- renderPlotly({
    req(volcano_data())
    df <- volcano_data()

    fc_thresh <- input$volcano_fc_threshold
    pval_thresh <- input$volcano_pval_threshold

    # Color palette
    if (input$volcano_color_mode == "viridis") {
      colors <- c("Up in Low Input" = "#440154", "Down in Low Input" = "#21918c", "Not Significant" = "#bbbbbb")
    } else {
      colors <- c("Up in Low Input" = "#e74c3c", "Down in Low Input" = "#3498db", "Not Significant" = "#bbbbbb")
    }

    p <- ggplot(df, aes(x = Log2FC, y = NegLog10P, color = Significant, text = Lipid)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = pval_thresh, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = colors) +
      labs(
        x = "Log2 Fold Change (Low Input / Control)",
        y = "-Log10(p-value)",
        color = "Significance"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")

    # Add labels for top significant lipids
    if (input$volcano_show_labels) {
      top_sig <- df %>%
        filter(Significant != "Not Significant") %>%
        arrange(desc(abs(Log2FC) * NegLog10P)) %>%
        head(input$volcano_top_labels)

      if (nrow(top_sig) > 0) {
        p <- p + geom_text(data = top_sig, aes(label = Lipid),
                           size = 3, vjust = -0.5, hjust = 0.5, check_overlap = TRUE)
      }
    }

    ggplotly(p, tooltip = c("text", "x", "y"))
  })

  output$volcano_summary <- renderUI({
    req(volcano_data())
    df <- volcano_data()

    n_up <- sum(df$Significant == "Up in Low Input")
    n_down <- sum(df$Significant == "Down in Low Input")
    n_total <- nrow(df)

    tags$div(
      style = "margin-top: 15px; padding: 10px; background: #f8f9fa; border-radius: 5px;",
      tags$strong("Summary: "),
      tags$span(class = "stat-badge significant", paste0(n_up, " Up")),
      tags$span(class = "stat-badge not-significant", paste0(n_down, " Down")),
      tags$span(style = "margin-left: 10px; color: #666;",
                paste0("(", n_total, " total lipids tested)"))
    )
  })

  output$volcano_table <- renderDT({
    req(volcano_data())
    df <- volcano_data() %>%
      arrange(PValue) %>%
      mutate(
        Log2FC = round(Log2FC, 3),
        PValue = signif(PValue, 3),
        NegLog10P = round(NegLog10P, 2),
        Mean_Control = round(Mean_Control, 2),
        Mean_LowInput = round(Mean_LowInput, 2)
      ) %>%
      select(Lipid, Class, Log2FC, PValue, NegLog10P, Mean_Control, Mean_LowInput, Significant)

    datatable(df, options = list(pageLength = 15, scrollX = TRUE), filter = "top")
  })

  # Enrichment Plot
  output$enrichment_plot <- renderPlot({
    req(volcano_data())
    df <- volcano_data()

    sig_lipids <- df %>% filter(Significant != "Not Significant")
    all_lipids <- df

    if (nrow(sig_lipids) < 5) {
      return(NULL)
    }

    # Count lipids per class in significant and all
    sig_counts <- table(sig_lipids$Class)
    all_counts <- table(all_lipids$Class)

    # Calculate enrichment (observed / expected)
    classes <- intersect(names(sig_counts), names(all_counts))

    enrichment_df <- data.frame(
      Class = classes,
      Significant = as.numeric(sig_counts[classes]),
      Total = as.numeric(all_counts[classes])
    ) %>%
      mutate(
        Expected = Total * (sum(Significant) / sum(Total)),
        Enrichment = log2(Significant / Expected),
        Direction = ifelse(Enrichment > 0, "Enriched", "Depleted")
      ) %>%
      filter(is.finite(Enrichment)) %>%
      arrange(desc(Enrichment))

    ggplot(enrichment_df, aes(x = reorder(Class, Enrichment), y = Enrichment, fill = Direction)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      coord_flip() +
      scale_fill_manual(values = c("Enriched" = "#27ae60", "Depleted" = "#e74c3c")) +
      labs(
        x = "Lipid Class",
        y = "Log2 Enrichment (Observed / Expected)",
        title = "Lipid Class Enrichment in Significant Lipids"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  })

  # Download handlers for volcano
  output$download_volcano_plot <- downloadHandler(
    filename = function() { "volcano_plot.png" },
    content = function(file) {
      df <- volcano_data()
      fc_thresh <- input$volcano_fc_threshold
      pval_thresh <- input$volcano_pval_threshold

      if (input$volcano_color_mode == "viridis") {
        colors <- c("Up in Low Input" = "#440154", "Down in Low Input" = "#21918c", "Not Significant" = "#bbbbbb")
      } else {
        colors <- c("Up in Low Input" = "#e74c3c", "Down in Low Input" = "#3498db", "Not Significant" = "#bbbbbb")
      }

      p <- ggplot(df, aes(x = Log2FC, y = NegLog10P, color = Significant)) +
        geom_point(alpha = 0.7, size = 2) +
        geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", color = "gray50") +
        geom_hline(yintercept = pval_thresh, linetype = "dashed", color = "gray50") +
        scale_color_manual(values = colors) +
        labs(x = "Log2 Fold Change (Low Input / Control)", y = "-Log10(p-value)") +
        theme_minimal(base_size = 14)

      ggsave(file, p, width = 10, height = 8, dpi = 300)
    }
  )

  output$download_volcano_data <- downloadHandler(
    filename = function() { "differential_lipids.csv" },
    content = function(file) {
      write.csv(volcano_data(), file, row.names = FALSE)
    }
  )

  ##======================================================================
  ##  NETWORK VISUALIZATION
  ##======================================================================

  network_data <- reactive({
    req(input$network_dataset, input$network_data_type, input$network_lipid_selection)

    df <- get_dataset(input$network_dataset, input$network_data_type)
    mat <- as.matrix(df[, -1])
    mat <- apply(mat, 2, as.numeric)

    # Select lipids
    if (input$network_lipid_selection == "top_var") {
      # Top variable lipids by MAD
      lipid_mad <- apply(mat, 2, function(x) mad(x, na.rm = TRUE))
      top_lipids <- names(sort(lipid_mad, decreasing = TRUE))[1:min(input$network_top_n, length(lipid_mad))]
      mat <- mat[, top_lipids, drop = FALSE]
    } else {
      # By lipid class
      lipid_classes <- sapply(colnames(mat), get_lipid_class)
      selected_lipids <- names(lipid_classes)[lipid_classes == input$network_class_filter]
      mat <- mat[, intersect(selected_lipids, colnames(mat)), drop = FALSE]
    }

    if (ncol(mat) < 5) {
      return(NULL)
    }

    # Calculate correlation matrix
    cor_mat <- cor(mat, use = "pairwise.complete.obs", method = input$network_method)

    # Create edge list
    threshold <- input$network_cor_threshold
    edges <- data.frame()

    for (i in 1:(ncol(cor_mat) - 1)) {
      for (j in (i + 1):ncol(cor_mat)) {
        cor_val <- cor_mat[i, j]
        if (!is.na(cor_val)) {
          if (abs(cor_val) >= threshold) {
            if (input$network_show_negative || cor_val > 0) {
              edges <- rbind(edges, data.frame(
                from = colnames(cor_mat)[i],
                to = colnames(cor_mat)[j],
                correlation = cor_val,
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
    }

    if (nrow(edges) == 0) {
      return(NULL)
    }

    # Create nodes
    all_nodes <- unique(c(edges$from, edges$to))
    node_degrees <- table(c(edges$from, edges$to))

    nodes <- data.frame(
      id = all_nodes,
      label = substr(all_nodes, 1, 15),  # Truncate long names
      title = all_nodes,  # Full name on hover
      group = sapply(all_nodes, get_lipid_class),
      value = as.numeric(node_degrees[all_nodes]),  # Node size by degree
      stringsAsFactors = FALSE
    )

    # Edge styling
    edges$color <- ifelse(edges$correlation > 0, "#27ae60", "#e74c3c")
    edges$width <- abs(edges$correlation) * 3
    edges$title <- paste0("r = ", round(edges$correlation, 3))

    list(nodes = nodes, edges = edges, cor_mat = cor_mat)
  })

  output$network_plot <- renderVisNetwork({
    data <- network_data()

    if (is.null(data)) {
      return(visNetwork(data.frame(id = 1, label = "Not enough data"), data.frame()) %>%
               visNodes(color = "#ccc"))
    }

    visNetwork(data$nodes, data$edges) %>%
      visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
        selectedBy = "group"
      ) %>%
      visGroups(groupname = "PC", color = "#00441B") %>%
      visGroups(groupname = "PE", color = "#41AB5D") %>%
      visGroups(groupname = "PG", color = "#78C679") %>%
      visGroups(groupname = "PA", color = "#1B7837") %>%
      visGroups(groupname = "PS", color = "#C2E699") %>%
      visGroups(groupname = "DG", color = "#54278F") %>%
      visGroups(groupname = "TG", color = "#ED804A") %>%
      visGroups(groupname = "MGDG", color = "#FBB4D9") %>%
      visGroups(groupname = "DGDG", color = "#F768A1") %>%
      visGroups(groupname = "SQDG", color = "#9D4D6C") %>%
      visGroups(groupname = "MG", color = "#8941ED") %>%
      visGroups(groupname = "Other", color = "#999999") %>%
      visPhysics(
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(gravitationalConstant = -50)
      ) %>%
      visInteraction(navigationButtons = TRUE, zoomView = TRUE) %>%
      visLegend(useGroups = TRUE, position = "right")
  })

  output$network_stats <- renderPrint({
    data <- network_data()

    if (is.null(data)) {
      cat("Not enough data to build network.\nTry lowering the correlation threshold.")
      return()
    }

    n_nodes <- nrow(data$nodes)
    n_edges <- nrow(data$edges)
    avg_degree <- mean(data$nodes$value)
    density <- (2 * n_edges) / (n_nodes * (n_nodes - 1))

    cat("Network Statistics\n")
    cat("==================\n")
    cat(sprintf("Nodes (Lipids): %d\n", n_nodes))
    cat(sprintf("Edges (Correlations): %d\n", n_edges))
    cat(sprintf("Average Degree: %.2f\n", avg_degree))
    cat(sprintf("Network Density: %.4f\n", density))
    cat(sprintf("Positive Edges: %d\n", sum(data$edges$correlation > 0)))
    cat(sprintf("Negative Edges: %d\n", sum(data$edges$correlation < 0)))
  })

  output$network_hubs <- renderDT({
    data <- network_data()

    if (is.null(data)) {
      return(datatable(data.frame(Message = "No network data available")))
    }

    hub_df <- data$nodes %>%
      arrange(desc(value)) %>%
      head(15) %>%
      mutate(
        Lipid = id,
        Class = group,
        Connections = value
      ) %>%
      select(Lipid, Class, Connections)

    datatable(hub_df, options = list(pageLength = 10, dom = 't'))
  })

  output$download_network_edges <- downloadHandler(
    filename = function() {
      paste0("network_edges_", input$network_dataset, ".csv")
    },
    content = function(file) {
      data <- network_data()
      if (!is.null(data)) {
        write.csv(data$edges[, c("from", "to", "correlation")], file, row.names = FALSE)
      }
    }
  )

  ##======================================================================
  ##  COMBINED DATA DOWNLOAD
  ##======================================================================
  output$download_combined_data <- downloadHandler(
    filename = function() { "combined_lipid_data.csv" },
    content = function(file) {
      df_control <- dataset_list$control
      df_lowinput <- dataset_list$lowinput

      all_lipids <- union(colnames(df_control)[-1], colnames(df_lowinput)[-1])
      ctrl_only <- setdiff(colnames(df_control)[-1], colnames(df_lowinput)[-1])
      low_only <- setdiff(colnames(df_lowinput)[-1], colnames(df_control)[-1])
      shared <- intersect(colnames(df_control)[-1], colnames(df_lowinput)[-1])

      combined_df <- data.frame(
        Lipid = all_lipids,
        Class = sapply(all_lipids, get_lipid_class),
        In_Control = all_lipids %in% colnames(df_control),
        In_LowInput = all_lipids %in% colnames(df_lowinput),
        Status = case_when(
          all_lipids %in% shared ~ "Shared",
          all_lipids %in% ctrl_only ~ "Control Only",
          all_lipids %in% low_only ~ "Low Input Only"
        )
      )

      write.csv(combined_df, file, row.names = FALSE)
    }
  )

}

shinyApp(ui, server)
