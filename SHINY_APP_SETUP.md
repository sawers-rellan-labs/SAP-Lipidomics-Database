# SAP Lipidomics Shiny App Setup

Run the app from the repository root with:

```r
shiny::runApp()
```

The app will start even if some files are missing, and it will show a setup banner listing what still needs to be added.

## Preferred file layout

Put these files inside the repository:

```text
app.R
data/
  raw/
    Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv
    Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv
  metadata/
    final_lipid_classes.csv
    gene_annotation.txt
  gwas/
    all_annotations_control.rds
    all_annotations_lowinput.rds
    all_annotations_control_PC.rds
    all_annotations_lowinput_PC.rds
www/
  manhattans/
    control/
      <one JPG per control lipid GWAS trait>
    lowinput/
      <one JPG per lowinput lipid GWAS trait>
    pca_control/
      <one JPG per control PCA GWAS trait>
    pca_lowinput/
      <one JPG per lowinput PCA GWAS trait>
```

## Required data files

These are required for the full app:

1. `Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv`
2. `Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv`
3. `final_lipid_classes.csv`
4. `gene_annotation.txt`
5. `all_annotations_control.rds`
6. `all_annotations_lowinput.rds`
7. `all_annotations_control_PC.rds`
8. `all_annotations_lowinput_PC.rds`

## Manhattan plot images

These are optional for startup, but required if you want the GWAS image panel to show the Manhattan plots instead of a placeholder message.

1. `www/manhattans/control/*.jpg`
2. `www/manhattans/lowinput/*.jpg`
3. `www/manhattans/pca_control/*.jpg`
4. `www/manhattans/pca_lowinput/*.jpg`

The JPG filename should match the GWAS trait name exactly, for example:

```text
www/manhattans/control/TG(53:6).jpg
www/manhattans/pca_lowinput/PC1_glycerolipid_li.jpg
```

## Alternate filenames the app also accepts

If you already have the older SoLD-style names, the app will also detect:

1. `data/Control_all_lipids_final.csv`
2. `data/Lowinput_all_lipids_final.csv`
3. `data/lipid_class.csv`
4. `data/all_annotations_control.rds`
5. `data/all_annotations_lowinput.rds`
6. `data/all_annotations_control_PC.rds`
7. `data/all_annotations_lowinput_PC.rds`

## R packages used by the app

Install these in R if needed:

```r
install.packages(c(
  "shiny", "bslib", "DT", "ggplot2", "dplyr", "vroom", "plotly",
  "Rtsne", "umap", "heatmaply", "RColorBrewer", "tibble",
  "VennDiagram", "ggnewscale", "stringr", "shinycssloaders",
  "visNetwork", "igraph", "viridisLite"
))
```
